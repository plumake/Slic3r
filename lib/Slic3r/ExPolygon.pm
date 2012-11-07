package Slic3r::ExPolygon;
use strict;
use warnings;

# an ExPolygon is a polygon with holes

use Boost::Geometry::Utils;
use Math::Geometry::Triangle;
use Slic3r::Geometry qw(X Y A B point_in_polygon same_line line_length distance_between_points);
use Slic3r::Geometry::Clipper qw(union_ex JT_MITER);
use List::Util qw(first);

# the constructor accepts an array of polygons 
# or a Math::Clipper ExPolygon (hashref)
sub new {
    my $class = shift;
    my $self;
    if (@_ == 1 && ref $_[0] eq 'HASH') {
        $self = [
            Slic3r::Polygon->new($_[0]{outer}),
            map Slic3r::Polygon->new($_), @{$_[0]{holes}},
        ];
    } else {
        $self = [ map Slic3r::Polygon->new($_), @_ ];
    }
    bless $self, $class;
    $self;
}

sub clone {
    my $self = shift;
    return (ref $self)->new(map $_->clone, @$self);
}

sub contour {
    my $self = shift;
    return $self->[0];
}

sub holes {
    my $self = shift;
    return @$self[1..$#$self];
}

sub lines {
    my $self = shift;
    return map $_->lines, @$self;
}

sub clipper_expolygon {
    my $self = shift;
    return {
        outer => $self->contour,
        holes => [ $self->holes ],
    };
}

sub boost_polygon {
    my $self = shift;
    return Boost::Geometry::Utils::polygon(@$self);
}

sub offset {
    my $self = shift;
    return Slic3r::Geometry::Clipper::offset($self, @_);
}

sub offset_ex {
    my $self = shift;
    return Slic3r::Geometry::Clipper::offset_ex($self, @_);
}

sub safety_offset {
    my $self = shift;
    
    # we're offsetting contour and holes separately
    # because Clipper doesn't return polygons in the same order as 
    # we feed them to it
    
    return (ref $self)->new(
        $self->contour->safety_offset,
        @{ Slic3r::Geometry::Clipper::safety_offset([$self->holes]) },
    );
}

sub noncollapsing_offset_ex {
    my $self = shift;
    my ($distance, @params) = @_;
    
    return $self->offset_ex($distance + 1, @params);
}

sub encloses_point {
    my $self = shift;
    my ($point) = @_;
    return $self->contour->encloses_point($point)
        && (!grep($_->encloses_point($point), $self->holes)
            || grep($_->point_on_segment($point), $self->holes));
}

# A version of encloses_point for use when hole borders do not matter.
# Useful because point_on_segment is slow
sub encloses_point_quick {
    my $self = shift;
    my ($point) = @_;
    return $self->contour->encloses_point($point)
        && !grep($_->encloses_point($point), $self->holes);
}

sub encloses_line {
    my $self = shift;
    my ($line, $tolerance) = @_;
    my $clip = $self->clip_line($line);
    if (!defined $tolerance) {
        # optimization
        return @$clip == 1 && same_line($clip->[0], $line);
    } else {
        return @$clip == 1 && abs(line_length($clip->[0]) - $line->length) < $tolerance;
    }
}

sub point_on_segment {
    my $self = shift;
    my ($point) = @_;
    for (@$self) {
        my $line = $_->point_on_segment($point);
        return $line if $line;
    }
    return undef;
}

sub bounding_box {
    my $self = shift;
    return Slic3r::Geometry::bounding_box($self->contour);
}

sub bounding_box_polygon {
    my $self = shift;
    my @bb = $self->bounding_box;
    return Slic3r::Polygon->new([
        [ $bb[0], $bb[1] ],
        [ $bb[2], $bb[1] ],
        [ $bb[2], $bb[3] ],
        [ $bb[0], $bb[3] ],
    ]);
}

sub bounding_box_center {
    my $self = shift;
    return Slic3r::Geometry::bounding_box_center($self->contour);
}

sub clip_line {
    my $self = shift;
    my ($line) = @_;  # line must be a Slic3r::Line object
    
    return Boost::Geometry::Utils::polygon_linestring_intersection(
        $self->boost_polygon,
        $line->boost_linestring,
    );
}

sub simplify {
    my $self = shift;
    $_->simplify(@_) for @$self;
}

sub translate {
    my $self = shift;
    $_->translate(@_) for @$self;
}

sub rotate {
    my $self = shift;
    $_->rotate(@_) for @$self;
}

sub area {
    my $self = shift;
    my $area = $self->contour->area;
    $area -= $_->area for $self->holes;
    return $area;
}

# Medial axis approximation based on Voronoi diagram.
# Returns a list of polylines and/or polygons.
# Trims leaves/twigs shorter than ($width/2)*sqrt(2).
sub medial_axis {
    my $self = shift;
    my ($width) = @_;
    
    my @self_lines = map $_->lines, @$self;
    my $expolygon = $self->clone;
    my @points = ();
    foreach my $polygon (@$expolygon) {
        Slic3r::Geometry::polyline_remove_short_segments($polygon, $width / 2);
        
        # subdivide polygon segments so that we don't have anyone of them
        # being longer than $width / 2
        $polygon->subdivide($width/2);
    }
    
    my $tri = Math::Geometry::Triangle->new();
    $tri->doVoronoi(1);
    $tri->addPolygon($expolygon->contour);
    $tri->addHole($_) for $expolygon->holes;
    my (undef, $vtopo) = $tri->triangulate('BEP'); # options supress unused output

    # remove references to ray edges from all nodes
    foreach my $node (@{$vtopo->{nodes}}) {
        # vector [0,0] means it's not a ray, so keep it
        @{$node->{edges}} = grep !$_->{vector}->[0] && !$_->{vector}->[1], @{$node->{edges}};
    }
    
    # all nodes where more than two edges meet are branch nodes
    my @branch_start_nodes = grep @{$_->{edges}} > 2, @{$vtopo->{nodes}};
    # if no branch nodes, dealing with a line or loop
    # if a line, nodes at ends will have only one edge reference 
    if (@branch_start_nodes == 0) {
        @branch_start_nodes = grep @{$_->{edges}} == 1, @{$vtopo->{nodes}};
    }
    # otherwise, it's a loop - any node can be the start node
    if (@branch_start_nodes == 0) {
        push @branch_start_nodes, $vtopo->{nodes}->[0];
    }

    my @polyedges = ();
    my @end_edges = ();
    
    # walk the cross referenced nodes and edges to build up polyline-like node lists
    foreach my $start_node (@branch_start_nodes) {
        foreach my $start_edge (@{$start_node->{edges}}) {
            # don't go back over path already traveled
            next if first {$_ == $start_edge} @end_edges;
            my $this_node = $start_node;
            push @polyedges, [];
            push @{$polyedges[-1]}, $this_node;
            my $this_edge = $start_edge;
            #step along nodes: next node is the node on current edge that isn't this one
            while ($this_node = +(grep $_ != $this_node, @{$this_edge->{nodes}})[0]) {
                # always add the point - duplicate start and end lets us detect polygons
                push @{$polyedges[-1]}, $this_node;
                # stop at a branch node, and remember the edge so we don't backtrack
                if (@{$this_node->{edges}} > 2) {
                    push @end_edges, $this_edge;
                    last;
                }
                # stop at a dead end
                last if @{$this_node->{edges}} == 1;
                # step to next edge
                $this_edge = +(grep $_ != $this_edge, @{$this_node->{edges}})[0];
                # stop if we've looped around to start
                last if $this_edge == $start_edge;
            }
        }
    }

    # If the medial_axis had branches, you might want to filter out leaf
    # edges. Or you might not. Not always clear when a leaf is a leaf.
    
    # This trims off any short polyedges that start at a branch node and end 
    # at a leaf node, where all nodes are within the half the length of the 
    # diagonal of a cube that encloses a circle of radius $width -
    # the idea being to trim off the leaf edges that go into the two corners
    # of the end of a rectangular area, since by the time we've reached the 
    # branch node those leaves are attached to, we're at $width/2 distance 
    # from the end of the rectange - the end.
    
    my $cuberad = ($width / 2) * sqrt(2);
    for (my $i = $#polyedges; $i > -1; $i--) {
        if (
            ( (@{$polyedges[$i]->[0]->{edges}}  == 1 && @{$polyedges[$i]->[-1]->{edges}} > 2) ||
              (@{$polyedges[$i]->[-1]->{edges}} == 1 && @{$polyedges[$i]->[0]->{edges}}  > 2)   )
            &&  distance_between_points($polyedges[$i]->[0]->{point}, $polyedges[$i]->[-1]->{point}) <= $cuberad
            ) {
            # candidate for trimming, based on leaf end node, 
            # so now check other nodes between start and end to be sure
            my $remove = 1;
            my $at_branch = @{$polyedges[$i]->[0]->{edges}} > 2 ? $polyedges[$i]->[0] : $polyedges[$i]->[-1];
            foreach my $p (@{$polyedges[$i]}[1 .. @{$polyedges[$i]} - 2]) {
                if (distance_between_points($at_branch->{point}, $p->{point}) > $cuberad) {
                    $remove = 0;
                    last;
                }
            }
            if ($remove) { splice(@polyedges, $i, 1); }
        }
    }
    
    # What about a branch going out to the tip of a long, narrow, 
    # tapering feature? Should we trim the tip back relative to the 
    # narrow tip angle? or relative to the distance-to-wall on either 
    # side of branch?
    # Distance-to-wall can probably be approximated by looking up edges 
    # in corresponding Delaunay triangulation. This would be close to a
    # "minimally inscribed circle" calculation for each medial axis node,
    # and could be refined.
    
    # Now combine chains of edges into longer chains where their ends meet,
    # deciding which two chains to link at branch node sites.

    my @polylines = ();

    # sort by length, where array length is rough proxy for edge length sum
    @polyedges = sort {@{$a} <=> @{$b}} @polyedges;

    # link polyedges with common end points
    for (my $i = $#polyedges; $i > 0; $i--) {
        # polygons
        if ($polyedges[$i]->[0] == $polyedges[$i]->[@{$polyedges[$i]} - 1]) {
            push @polylines, splice(@polyedges, $i, 1);
            next;
        }
        # polylines
        else {
            # We take a longer polyline from the end of the list
            # and see if it links up with any of the shorter
            # polylines that come before it in the list.
            # If so, we splice the longer polyline out of the list
            # and add it's points to the shorter polyline.
            # Having that first splice choose the shorter next path
            # is meant to fill out local features while in the neighborhood
            # instead of always linking longest paths, which might
            # require revisting a lot of separate regions of local features
            # that were passed by, requiring more rapid traversals.
            # The paths may need further sorting after this linking to achieve
            # this. The point here is just to link up the path structure that 
            # will enable that.
            # ... sort of works - but there's probably a better approach
            
            my $this  = $polyedges[$i];

            for (my $j = 0; $j < $i ; $j++) {
                my $other = $polyedges[$j];
                # all the cases of ends matching up
                if ($this->[@{$this} - 1] == $other->[0]) {
                    shift @{$other};
                    @{$other} = (@{splice(@polyedges, $i, 1)}, @{$other});
                    last;
                } elsif ($this->[0] == $other->[@{$other} - 1]) {
                    shift @{$this};
                    @{$other} = (@{$other}, @{splice(@polyedges, $i, 1)});
                    last;
                } elsif ($this->[0] == $other->[0]) {
                    shift @{$this};
                    @{$other} = ((reverse @{$other}), @{splice(@polyedges, $i, 1)});
                    last;
                } elsif ($this->[@{$this} - 1] == $other->[@{$other} - 1]) {
                    pop @{$other};
                    @{$other} = (@{splice(@polyedges, $i ,1)}, (reverse @{$other}));
                    last;
                }
            }
        }
    }

    push @polylines, @polyedges;
    
    my @result = ();
    foreach my $polyline (@polylines) {
        next unless @$polyline >= 2;
        
        # now extract just the point coordinates from the nodes
        @$polyline = map $_->{point}, @$polyline;
                
        # cleanup
        #$polyline = Slic3r::Geometry::douglas_peucker($polyline, $width / 7);
        
        if (Slic3r::Geometry::same_point($polyline->[0], $polyline->[-1])) {
            next if @$polyline == 2;
            push @result, Slic3r::Polygon->new(@$polyline[0..$#$polyline-1]);
        } else {
            push @result, Slic3r::Polyline->new($polyline);
        }
    }
    
    return @result;
}

1;
