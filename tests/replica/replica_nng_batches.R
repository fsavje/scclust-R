batch_assign <- function(nng,
                         unassigned_method) {
  cl_label <- rep(as.integer(NA), ncol(nng))
  current_label <- 0
  assigned <- rep(FALSE, ncol(nng))
  invisible(lapply(1:ncol(nng), function(i) {
    if (!assigned[i] && any(nng[, i])) {
      if (!any(assigned[nng[, i]])) {
        assigned[nng[, i]] <<- TRUE
        cl_label[nng[, i]] <<- current_label
        current_label <<- current_label + 1
      } else if (unassigned_method == "any_neighbor") {
        cand <- which(nng[, i])
        cl_label[i] <<- cl_label[cand[assigned[cand]][1]]
      }
    }
  }))

  cl_label
}


replica_nng_clustering_batches <- function(distances,
                                           size_constraint,
                                           unassigned_method = "any_neighbor",
                                           radius = NULL,
                                           primary_data_points = NULL) {
  ensure_distances(distances)
  num_data_points <- length(distances)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  unassigned_method <- coerce_args(unassigned_method,
                                   c("ignore", "any_neighbor"))
  radius <- coerce_radius(radius)
  if (!is.null(primary_data_points)) {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }

  distances <- as.matrix(distances)
  nng <- get_simple_nng(distances,
                        size_constraint,
                        radius,
                        primary_data_points)

  cl_label <- batch_assign(nng, unassigned_method)

  make_scclust(as.integer(cl_label),
               length(unique(cl_label[!is.na(cl_label)])),
               attr(distances, "ids", exact = TRUE))
}
