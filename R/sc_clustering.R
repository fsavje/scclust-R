# ==============================================================================
# scclust for R -- R wrapper for the scclust library
# https://github.com/fsavje/scclust-R
#
# Copyright (C) 2016  Fredrik Savje -- http://fredriksavje.com
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


#' Size constrained clustering
#'
#' \code{sc_clustering} constructs near-optimal size constrained clusterings.
#' Subject to user-specified conditions on the minimum size and composition
#' of the clusters, the function derives a partition of a set of data points
#' so that the dissimilarity of points assigned to the same cluster is
#' minimized.
#'
#' \code{sc_clustering} implements an algorithm that first summaries the distance
#' information between data points in a sparse graph and then constructs the clusterings
#' based on the graph. This makes the function run fast while ensuring
#' near-optimal performance. It is possible to constrain the clustering so that each cluster contains at
#' least a certain number of points in total. If there are different types of points, each cluster can also be
#' constrained to contain a certain number of points of
#' each type. For example, in a sample with
#' "red" and "blue" data points, one can constrain the clusters so that each contains
#' at least 10 points in total and at least 5 "blue" points.
#'
#' The vertices in the graph that \code{sc_clustering} constructs represent data points. Arcs
#' in the graph are weighted by the distance between the data points.
#' The graph is the smallest graph such that the neighborhood of each vertex satsifies
#' the clustering constraints supplied by the user. By picking vertices ("seeds") with
#' non-overlapping neighborhoods and constructing clusters as supersets of the neighborhoods,
#' we ensure that the clustering will satisfy the constraints. While all methods of picking seeds yield near-optimal clusterings, we typically want as
#' many seeds as possible; this will lead to many clusters and, thus, smaller clusters. The \code{seed_method} option specifies how the seeds are picked.
#'
#' When \code{seed_method} is set to "lexical", seeds are picked in lexical order. This
#' allows us to find the seeds quickly. However, we might pick points
#' that are central in the sparse
#' graph with this method. When central points are seeds, they will exclude many other points from being seeds and, thus, lead to larger clusters The "exclusion_order" and "exclusion_updating" options
#' calculates, for each vertex, how many other vertices are excluded
#' if the vertex is picked. By picking vertices that exclude few vertices, we avoid
#' central points and, thereby, increase
#' the number of seeds. "exclusion_updating" updates the count after each picked seed
#' so that already excluded vertices are not counted twice; "exclusion_order" derives
#' the count once.
#'
#' Deriving the exclusion count is an expensive operation. The
#' "inwards_order" and "inwards_updating" options count the number of inwards-pointing
#' arcs in the graph, which approximates the exclusion count.  "inwards_updating" updates
#' the count after each picked seed, while "inwards_order" derives the count once. The "batches" option is identical to "lexical" but it derives the graph in batches.
#' This limits memory use to a constant value decided by \code{batch_size}. This can be
#' useful with large size constraints. The "batches" option is still experimental and
#' can currently only be used when there is no type constraints.
#'
#' The \code{seed_radius} option limits the maximum distance between adjencent vertices
#' in the sparse graph. If a distance in a point's neighborhood is greater than the radius,
#' the neighborhood will be deleted and the point cannot be a seed (it can, however, be in
#' other points' neighborhoods). As
#'
#' Some data points might not be in a neighborhood of a seed and will, therefore, not be
#' assigned to a cluster at this point. \code{primary_unassigned_method} specifies how
#' the unassigned points are assigned. With the "ignore" option, the points are left unassigned.
#' With the "any_neighbor", the function re-uses the sparse graph and assigns the unassigned
#' points to any neighborhood that is adjencent to the point in the graph. The "closest_assigned"
#' option assigns each point to the cluster that contain its closest assigned vertex, and the
#' "closest_seed" option assigns to the cluster that contain the closest seed. All these options
#' ensure near-optimal clusterings.
#'
#' Radius ..
#'
#' Occasionally, some data points in our clustering problem are allowed to be left unassigned.
#' For example, assume that we are using clustering for a prediction problem. We have a set of
#' training points with known values of some outcome of interest, and a set of test points
#' whose values we want to predict. Clustering the points, we can use the cluster mean of the
#' training point to predict the values of the test points in the cluster. To make this viable,
#' we need that each cluster contains at least one training point -- this is exactly what
#' \code{type_labels} and \code{type_constraints} do for us. We also need that each test point
#' is assigned to a cluster. We do, however, not need that each training point is assigned.
#' In this case, by specifying only the test points in \code{primary_data_points}, we ensure
#' that each of these points will be assigned to clusters. The points not specified in
#' \code{primary_data_points} (i.e., the training points in our case) will be assigned to
#' clustering only insofar that it is needed to satisfy the clustering constraints. These
#' points are called "secondary". This can lead to large improvements in the clustering
#' if the types of points are unevenly spaced.
#'
#' In most cases, we do not want to discard all unassigned secondary points -- some of them
#' will, occasionally, be close to a cluster and contain useful information. Similar to
#' \code{primary_unassigned_method}, we can use the \code{secondary_unassigned_method}
#' to specify how the leftover secondary points should be assigned. The three possible
#' options are "ignore", "closest_assigned" and "closest_seed".
#'
#' Radius...
#'
#'
#'
#' @param distances a \code{\link[distances]{distances}} object with distances
#'                  between the data points.
#'
#' @param size_constraint an integer with the required minimum cluster size.
#'
#' @param type_labels a vector containing the type of each data point. May be
#'                    \code{NULL} when \code{type_constraints} is \code{NULL}.
#'
#' @param type_constraints a named integer vector containing type-specific size constraints.
#'                         If \code{NULL}, only the overall constraint given by
#'                         \code{size_constraint} will be imposed on the clustering.
#'
#' @param seed_method a character scalar indicating how seeds should be found.
#'
#' @param primary_data_points a vector specifying primary data points, either by point indices
#'                            or with a logical vector of length equal to the number of points.
#'                            \code{NULL} indicates that all data points are "primary".
#'
#' @param primary_unassigned_method a character scalar indicating how unassigned (primary) points
#'                                  should be assigned to clusters.
#'
#' @param secondary_unassigned_method a character scalar indicating how unassigned secondary points
#'                                    should be assigned to clusters.
#'
#' @param seed_radius a positive numeric scalar restricting the maximum length of an edge in the
#'                    graph used to construct the clustering. \code{NULL} indicates
#'                    no restriction. This parameter (together with \code{primary_radius} and
#'                    \code{secondary_radius}) can be used to control the maximum distance
#'                    between points assigned to the same cluster (see below for details).
#'
#' @param primary_radius a positive numeric scalar, a character scalar or NULL restricting the
#'                       match distance for unassigned primary points. If numeric, the value
#'                       is used to restrict the distances. If character, it must be one of
#'                       "no_radius", "seed_radius" or "estimated_radius" (see below for
#'                       details). \code{NULL} indicates no radius restriction.
#'
#' @param secondary_radius a positive numeric scalar, a character scalar or NULL restricting the
#'                         match distance for unassigned primary points. If numeric, the value
#'                         is used to restrict the distances. If character, it must be one of
#'                         "no_radius", "seed_radius" or "estimated_radius" (see below for
#'                         details). \code{NULL} indicates no radius restriction.
#'
#' @param batch_size an integer scalar specifying batch size when \code{seed_method} is set to
#'                   "batches".
#'
#' @return Returns a \code{\link{scclust}} object with the derived clustering.
#'
#' @seealso \code{\link{hierarchical_clustering}} can be used to refine the clustering
#'          constructed by \code{sc_clustering}.
#'
#' @examples
#' # Make example data
#' my_data <- data.frame(id = 1:100000,
#'                       type = factor(rbinom(100000, 3, 0.3),
#'                                     labels = c("A", "B", "C", "D")),
#'                       x1 = rnorm(100000),
#'                       x2 = rnorm(100000),
#'                       x3 = rnorm(100000))
#'
#' # Construct distance metric
#' my_dist <- distances(my_data,
#'                      id_variable = "id",
#'                      dist_variables = c("x1", "x2", "x3"),
#'                      normalize = "mahalanobize")
#'
#' # Make clustering with at least 3 data points in each cluster
#' my_clustering <- sc_clustering(my_dist, 3)
#'
#' # Check so clustering satisfies constraints
#' check_scclust(my_clustering, 3)
#' # > TRUE
#'
#' # Get statistics about the clustering
#' get_scclust_stats(my_clustering, my_dist)
#' # > num_data_points        1.000000e+05
#' # > ...
#'
#' # Make clustering with at least one point of each type in each cluster
#' my_clustering <- sc_clustering(my_dist,
#'                                type_labels = my_data$type,
#'                                type_constraints = c("A" = 1, "B" = 1,
#'                                                     "C" = 1, "D" = 1))
#'
#' # Check so clustering satisfies constraints
#' check_scclust(my_clustering,
#'               type_labels = my_data$type,
#'               type_constraints = c("A" = 1, "B" = 1,
#'                                    "C" = 1, "D" = 1))
#' # > TRUE
#'
#' # Make clustering with at least 8 points in total of which at least
#' # one must be "A" and 2 "Bs" and five can be any type
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
                          seed_method = "exclusion_updating",
                          primary_data_points = NULL,
                          primary_unassigned_method = "closest_seed",
                          secondary_unassigned_method = "ignore",
                          seed_radius = NULL,
                          primary_radius = "seed_radius",
                          secondary_radius = NULL,
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
