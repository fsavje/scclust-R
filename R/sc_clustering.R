# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016-2017  Fredrik Savje -- http://fredriksavje.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/
# ==============================================================================


#' Size-constrained clustering
#'
#' \code{sc_clustering} constructs near-optimal size-constrained clusterings.
#' Subject to user-specified constraints on the size and composition of the
#' clusters, a clustering is constructed so that within-cluster pair-wise
#' distances are minimized. The function does not restrict the number of
#' clusters.
#'
#' \code{sc_clustering} constructs a clustering so to minimize within-cluster
#' dissimilarities while ensuring that constraints on the size and composition
#' of the clusters are satisfied. It is possible to impose an overall size
#' constraint so that each cluster must contain at least a certain number of
#' points in total. It is also possible to impose constraints on the composition
#' of the clusters so that each cluster must contain a certain number of points
#' of different types. For example, in a sample with "red" and "blue" data
#' points, one can constrain the clustering so that each cluster must contain
#' at least 10 points in total of which at least 3 must be "red" and at least
#' 2 must be "blue".
#'
#' The function implements an algorithm that first summaries the distances
#' between data points in a sparse graph and then constructs the clustering
#' based on the graph. This admits fast execution while ensuring near-optimal
#' performance. In particular, the maximum within-cluster distance is garantueed
#' to be at most four times the maximum distance in the optimal clustering.
#' The average performance is much closer to optimal than the worst case bound.
#'
#' In more detail, the clustering algorithm has four steps:
#'
#' \enumerate{
#'   \item Construct sparse graph encoding clustering constraints.
#'   \item Select a set of vertices in the graph to act as "seeds".
#'   \item Construct initial clusters from the seeds' neighborhoods in the graph.
#'   \item Assign remaining data points to clusters.
#' }
#'
#' Each data point is represented by a vertex in the sparse graph, and arcs are
#' weighted by the distance between the data points they connect. The graph is
#' constructed so that each vertex's neighborhood satisfies the clustering
#' constraints supplied by the user. For example, if the clusters must contain
#' at least five points, each closed neighborhood in the graph will contain five
#' vertices. \code{sc_clustering} constructs the graph that minimize the arc
#' weights subject to the neighborhood constraints; this ensures near-optimal
#' performance. The function selects "seeds" that have non-overlapping
#' neighborhood in the graph. By constructing clusters as supersets of the
#' neighborhoods, it ensures that the clustering will satisfy the constraints
#' imposed by the user.
#'
#' The \code{seed_method} option governs how the seeds are chosen. Any set of
#' seeds yields a near-optimal clustering, but, heuristically, performance can
#' be improved by picking the seeds more carefully. In most cases, smaller
#' clusters are desirable since they tend to minimize within-cluster distances.
#' As the number of data points is fixed, we minimize the cluster size by
#' maximizing the number of clusters (and, thus, the number of seeds). When
#' \code{seed_method} is set to "lexical", seeds are chosen in lexical order.
#' This admits a fast solution, but the function might pick points that are
#' central in the graph. Central points tend to exclude many other points from
#' being seeds and, thus, lead to larger clusters.
#'
#' The "exclusion_order" and "exclusion_updating" options calculate for each
#' vertex how many other vertices are excluded if the vertex is picked as a seed.
#' By picking seeds that exclude few other vertices, the function avoids
#' central points and increases the number of seeds. "exclusion_updating"
#' updates the count after each picked seed so that already excluded vertices
#' are not counted twice; "exclusion_order" derives the count once. The former
#' option is, thus, better-performing but slower.
#'
#' Deriving the exclusion count when using the "exclusion_order" and "exclusion_updating"
#' options is an expensive operation, and it might not be feasible to do so in
#' large samples. This is particularly problematic when the data points have
#' equidistant nearest neighbors, which tends to happen when the dataset contains
#' only discrete variables. It also happens when the dataset contains many
#' identical data points. The exclusion count operation may in these cases take
#' several orders of magnitude longer to run than the rest of the operations
#' combined. To ensure sane behavior for the default options, the exclusion
#' count is not used by default. If the dataset contains at least one continuous
#' variable, it is generally safe to call \code{sc_clustering} with either
#' "exclusion_order" or "exclusion_updating", and this will often improve performance
#' over the default.
#'
#' The "inwards_order" and "inwards_updating" options count the number of
#' inwards-pointing arcs in the graph, which approximates the exclusion count.
#' "inwards_updating" updates the count after each picked seed, while
#' "inwards_order" derives the count once. The inwards counting options work
#' well with nearly all types of data, and "inwards_updating" is the default.
#'
#' The "batches" option is identical to "lexical" but it derives the graph in
#' batches. This limits the use of memory to a value proportional to
#' \code{batch_size} irrespectively of the size of the dataset. This can be
#' useful when imposing large size constraints, which consume a lot of memory
#' in large datasets. The "batches" option is still experimental and can currently
#' only be used when one does not impose type constraints.
#'
#' Once the function has selected seeds so that no additional points can be
#' selected without creating overlap in the graph, it constructs the initial
#' clusters. A unique cluster label is assigned to each seed, and all points
#' in the seed's neighborhood is assigned the same label. Some
#' data points might not be in a neighborhood of a seed and will, therefore,
#' not be assigned to a cluster after the initial clusters have been formed.
#' \code{primary_unassigned_method} specifies how the unassigned points are
#' assigned. With the "ignore" option, the points are left unassigned. With the
#' "any_neighbor", the function re-uses the sparse graph and assigns the
#' unassigned points to any adjacent cluster in the graph. The "closest_assigned"
#' option assigns the points to the clusters that contain their closest assigned
#' vertex, and the "closest_seed" option assigns to the cluster that contain
#' the closest seed. All these options ensure that the clustering is near-optimal.
#'
#' Occasionally, some data points are allowed to be
#' left unassigned.  Consider the following prediction problem as an example. We
#' have a set of data points with a known outcome value ("training points") and
#' another set of points for which the outcome is unknown ("prediction points").
#' We want to use the training points to predict the outcome for the prediction
#' points. By clustering the data, we can use the cluster mean of the training
#' points to predict the values of the prediction points in the same cluster. To
#' make this viable, we need to ensure that each cluster contains at least one
#' training point (the cluster mean would otherwise be undefined). We can impose
#' this constraint using the \code{type_labels} and \code{type_constraints}
#' options. We also need to make sure that each prediction point is assigned to
#' a cluster. We do, however, not need that all training points are assigned a
#' cluster; some training points might not provide useful information (e.g.,
#' if they are very dissimilar to all prediction points). In this case, by
#' specifying only the prediction points in \code{primary_data_points}, we
#' ensure that all those points are assigned to clusters. The points not
#' specified in \code{primary_data_points} (i.e., the training points) will
#' be assigned to clustering only insofar that it is needed to satisfy the
#' clustering constraints. This can lead to large improvements in the clustering
#' if the types of points are unevenly distributed in the metric space. Points
#' not specified in \code{primary_data_points} are called "secondary".
#'
#' Generally, one does not want to discard all unassigned secondary points --
#' some of them will, occasionally, be close to a cluster and contain useful
#' information. Similar to \code{primary_unassigned_method}, we can use the
#' \code{secondary_unassigned_method} to specify how the leftover secondary
#' points should be assigned. The three possible options are "ignore",
#' "closest_assigned" and "closest_seed". In nearly all cases, it is beneficial
#' to impose a radius constraint when assigning secondary points (see below for
#' details).
#'
#' \code{sc_clustering} tries to minimize within-cluster distances subject to
#' that all (primary) data points are assigned to clusters. In some cases, it
#' might be beneficial to leave some points unassigned to avoid clusters with
#' very dissimilar points. The function has options that can be used to
#' indirectly constrain the maximum within-cluster distance. Specifically, by
#' restricting the maximum distance in the four steps of the function, we can
#' bound the maximum within-cluster distance. The bound is, however, a blunt
#' tool. A too small bound might lead to that only a few data points are
#' assigned to clusters. If a bound is needed, it is recommended to use a very
#' liberal bound. This will avoid the very dissimilar clusters, but let the
#' function optimize the remaining cluster assignments.
#'
#' The \code{seed_radius} option limits the maximum distance between adjacent
#' vertices in the sparse graph. If the distance to a neighbor in a point's
#' neighborhood is greater than \code{seed_radius}, the neighborhood will be
#' deleted and the point is not allowed be a seed (it can, however, still be
#' in other points' neighborhoods). The \code{primary_radius} option limits
#' the maximum distance when a primary point is assigned to in the fourth step.
#' When \code{primary_radius} is set to a positive numeric value, this will be
#' used to as the radius restriction. "no_radius" and \code{NULL} indicate no
#' restriction. "seed_radius" indicates that the same restriction as
#' \code{seed_radius} should be used, and "estimated_radius" sets the
#' restriction to the estimated average distance between the seeds and their
#' neighbors. When \code{primary_unassigned_method} is set to "any_neighbor",
#' \code{primary_radius} must be set to "seed_radius".
#'
#' The way the radius constraints restrict the maximum within-cluster distance
#' depends on the \code{primary_unassigned_method} option. When
#' \code{primary_unassigned_method} is "ignore", the maximum distance is bounded
#' by 2 * \code{seed_radius}. When \code{primary_unassigned_method} is
#' "any_neighbor", it is bounded by 4 * \code{seed_radius}. When
#' \code{primary_unassigned_method} is "closest_assigned", it is bounded by
#' 2 * \code{seed_radius} + 2 * \code{primary_radius}. When
#' \code{primary_unassigned_method} is "closest_seed", it is bounded by the
#' maximum of 2 * \code{seed_radius} and 2 * \code{primary_radius}.
#'
#' The \code{secondary_radius} option restricts how secondary points are
#' assigned, but is otherwise identical to the \code{primary_radius} option.
#'
#' @param distances
#'    a \code{\link[distances]{distances}} object with distances between the
#'    data points.
#' @param size_constraint
#'    an integer with the required minimum cluster size.
#' @param type_labels
#'    a vector containing the type of each data point. May be \code{NULL} when
#'    \code{type_constraints} is \code{NULL}.
#' @param type_constraints
#'    a named integer vector containing type-specific size constraints. If
#'    \code{NULL}, only the overall constraint given by \code{size_constraint}
#'    will be imposed on the clustering.
#' @param seed_method
#'    a character scalar indicating how seeds should be selected.
#' @param primary_data_points
#'    a vector specifying primary data points, either with point indices or with
#'    a logical vector of length equal to the number of points. \code{NULL}
#'    indicates that all data points are "primary".
#' @param primary_unassigned_method
#'    a character scalar indicating how unassigned (primary) points should be
#'    assigned to clusters.
#' @param secondary_unassigned_method
#'    a character scalar indicating how unassigned secondary points should be
#'    assigned to clusters.
#' @param seed_radius
#'    a positive numeric scalar restricting the maximum length of an edge in the
#'    graph used to construct the clustering. \code{NULL} indicates no restriction.
#'    This parameter (together with \code{primary_radius} and \code{secondary_radius})
#'    can be used to control the maximum distance between points assigned to the
#'    same cluster (see below for details).
#' @param primary_radius
#'    a positive numeric scalar, a character scalar or NULL restricting the match
#'    distance for unassigned primary points. If numeric, the value is used to
#'    restrict the distances. If character, it must be one of "no_radius",
#'    "seed_radius" or "estimated_radius" (see below for details). \code{NULL}
#'    indicates no radius restriction.
#' @param secondary_radius
#'    a positive numeric scalar, a character scalar or NULL restricting the match
#'    distance for unassigned secondary points. If numeric, the value is used to
#'    restrict the distances. If character, it must be one of "no_radius",
#'    "seed_radius" or "estimated_radius" (see below for details). \code{NULL}
#'    indicates no radius restriction.
#' @param batch_size
#'    an integer scalar specifying batch size when \code{seed_method} is set to
#'    "batches".
#'
#' @return
#'    Returns a \code{\link{scclust}} object with the derived clustering.
#'
#' @seealso
#'    \code{\link{hierarchical_clustering}} can be used to refine the clustering
#'    constructed by \code{sc_clustering}.
#'
#' @references
#' Higgins, Michael J., Fredrik Sävje and Jasjeet S. Sekhon (2016),
#' \sQuote{Improving massive experiments with threshold blocking},
#' \emph{Proceedings of the National Academy of Sciences}, \bold{113:27}, 7369--7376.
#'
#' Sävje, Fredrik and Michael J. Higgins and Jasjeet S. Sekhon (2017),
#' \sQuote{Generalized Full Matching}, arXiv 1703.03882.
#' \url{https://arxiv.org/abs/1703.03882}
#'
#' @examples
#' # Make example data
#' my_data <- data.frame(id = 1:50000,
#'                       type = factor(rbinom(50000, 3, 0.3),
#'                                     labels = c("A", "B", "C", "D")),
#'                       x1 = rnorm(50000),
#'                       x2 = rnorm(50000),
#'                       x3 = rnorm(50000))
#'
#' # Construct distance metric
#' my_dist <- distances(my_data,
#'                      id_variable = "id",
#'                      dist_variables = c("x1", "x2", "x3"))
#'
#' # Make clustering with at least 3 data points in each cluster
#' my_clustering <- sc_clustering(my_dist, 3)
#'
#' # Check so clustering satisfies constraints
#' check_clustering(my_clustering, 3)
#' # > TRUE
#'
#' # Get statistics about the clustering
#' get_clustering_stats(my_dist, my_clustering)
#' # > num_data_points        5.000000e+04
#' # > ...
#'
#' # Make clustering with at least one point of each type in each cluster
#' my_clustering <- sc_clustering(my_dist,
#'                                type_labels = my_data$type,
#'                                type_constraints = c("A" = 1, "B" = 1,
#'                                                     "C" = 1, "D" = 1))
#'
#' # Check so clustering satisfies constraints
#' check_clustering(my_clustering,
#'                  type_labels = my_data$type,
#'                  type_constraints = c("A" = 1, "B" = 1,
#'                                       "C" = 1, "D" = 1))
#' # > TRUE
#'
#' # Make clustering with at least 8 points in total of which at least
#' # one must be "A", two must be "B" and five can be any type
#' my_clustering <- sc_clustering(my_dist,
#'                                size_constraint = 8,
#'                                type_labels = my_data$type,
#'                                type_constraints = c("A" = 1, "B" = 2))
#'
#' @keywords cluster
#' @export
sc_clustering <- function(distances,
                          size_constraint = NULL,
                          type_labels = NULL,
                          type_constraints = NULL,
                          seed_method = "inwards_updating",
                          primary_data_points = NULL,
                          primary_unassigned_method = "closest_seed",
                          secondary_unassigned_method = "ignore",
                          seed_radius = NULL,
                          primary_radius = "seed_radius",
                          secondary_radius = "estimated_radius",
                          batch_size = 100L) {
  ensure_distances(distances)
  num_data_points <- length(distances)

  if (is.null(type_constraints)) {
    type_labels <- NULL
    size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  } else {
    type_labels <- coerce_type_labels(type_labels, num_data_points)
    type_constraints <- coerce_type_constraints(type_constraints)
    type_constraints <- make_type_size_constraints(type_constraints,
                                                   type_labels)
    size_constraint <- coerce_total_size_constraint(size_constraint,
                                                    type_constraints,
                                                    num_data_points)
  }

  seed_method <- coerce_args(seed_method, all_seed_methods)

  primary_data_points <- coerce_data_point_indices(primary_data_points, num_data_points)
  if (is.logical(primary_data_points)) {
    primary_data_points <- which(primary_data_points) - 1L
  } else if (is.integer(primary_data_points)) {
    primary_data_points <- primary_data_points - 1L
  }

  primary_unassigned_method <- coerce_args(primary_unassigned_method,
                                           c("ignore",
                                             "any_neighbor",
                                             "closest_assigned",
                                             "closest_seed"))
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed"))

  seed_radius <- coerce_radius(seed_radius, TRUE)
  primary_radius <- coerce_radius(primary_radius, FALSE)
  secondary_radius <- coerce_radius(secondary_radius, FALSE)

  if (!is.null(batch_size)) {
    batch_size <- coerce_counts(batch_size, 1L)
  }

  clustering <- .Call(Rscc_sc_clustering,
                      distances,
                      size_constraint,
                      unclass(type_labels),
                      type_constraints,
                      seed_method,
                      primary_data_points,
                      primary_unassigned_method,
                      secondary_unassigned_method,
                      seed_radius,
                      primary_radius,
                      secondary_radius,
                      batch_size)

  make_scclust(clustering$cluster_labels,
               clustering$cluster_count,
               attr(distances, "ids", exact = TRUE))
}
