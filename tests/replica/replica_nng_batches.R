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
      } else if (unassigned_method == "by_nng") {
        cand <- which(nng[, i])
        cl_label[i] <<- cl_label[cand[assigned[cand]][1]]
      }
    }
  }))

  cl_label
}


replica_nng_clustering_batches <- function(distance_object,
                                           size_constraint,
                                           unassigned_method = "by_nng",
                                           radius = NULL,
                                           primary_data_points = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  unassigned_method <- coerce_args(unassigned_method,
                                   c("ignore", "by_nng"))
  radius <- coerce_radius(radius)
  if (!is.null(primary_data_points)) {
    ensure_indicators(primary_data_points, num_data_points, TRUE)
  }

  distances <- as.matrix(distance_object)
  nng <- get_simple_nng(distances,
                        size_constraint,
                        radius,
                        primary_data_points)

  cl_label <- batch_assign(nng, unassigned_method)

  make_Rscc_clustering(as.integer(cl_label),
                       length(unique(cl_label[!is.na(cl_label)])),
                       attr(distance_object, "ids", exact = TRUE))
}
