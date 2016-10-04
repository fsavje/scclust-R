batch_assign <- function(nng,
                         main_unassigned_method) {
  cl_label <- rep(as.integer(NA), ncol(nng))
  current_label <- 0
  assigned <- rep(FALSE, ncol(nng))
  invisible(lapply(1:ncol(nng), function(i) {
    if (!assigned[i] && any(nng[, i])) {
      if (!any(assigned[nng[, i]])) {
        assigned[nng[, i]] <<- TRUE
        cl_label[nng[, i]] <<- current_label
        current_label <<- current_label + 1
      } else if (main_unassigned_method == "by_nng") {
        cand <- which(nng[, i])
        cl_label[i] <<- cl_label[cand[assigned[cand]][1]]
      }
    }
  }))

  cl_label
}


replica_nng_clustering_batches <- function(distance_object,
                                           size_constraint,
                                           main_unassigned_method = "by_nng",
                                           main_radius = NULL,
                                           main_data_points = NULL) {
  check_Rscc_distances(distance_object)
  num_data_points <- get_num_data_points(distance_object)
  size_constraint <- get_size_constraint(size_constraint)
  stopifnot(size_constraint <= num_data_points)
  main_unassigned_method <- match.arg(main_unassigned_method, c("ignore", "by_nng"))
  main_radius <- get_radius(main_radius)
  check_main_data_points(main_data_points, num_data_points)

  distances <- as.matrix(distance_object)
  nng <- get_simple_nng(distances,
                        size_constraint,
                        main_radius,
                        main_data_points)

  cl_label <- batch_assign(nng, main_unassigned_method)

  make_Rscc_clustering(as.integer(cl_label),
                       length(unique(cl_label[!is.na(cl_label)])),
                       attr(distance_object, "ids", exact = TRUE))
}
