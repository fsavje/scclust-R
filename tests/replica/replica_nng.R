assign_by_seed <- function(seeds, nng) {
  cl_label <- rep(as.integer(NA), ncol(nng))
  current_label <- 0

  invisible(lapply(seeds, function(i) {
    cl_label[nng[, i]] <<- current_label
    current_label <<- current_label + 1
  }))

  cl_label
}


match_n_assign <- function(cl_label,
                           match,
                           unassigned,
                           radius,
                           distances) {
  if (!is.null(radius)) {
    unassigned[unassigned] <- (apply(distances[unassigned, match, drop = FALSE], 1, sort)[1, ] <= radius)
  }
  if (sum(unassigned) > 0) {
    cl_label[unassigned] <- cl_label[match[apply(distances[unassigned, match, drop = FALSE], 1, order)[1, ]]]
  }
  cl_label
}


by_nng_or_closest_assigned <- function(main_unassigned_method,
                                       main_unassigned,
                                       nng,
                                       assigned,
                                       main_radius,
                                       cl_label,
                                       distances) {
  if (main_unassigned_method == "by_nng") {
    main_unassigned <- apply(nng, 2, any) & main_unassigned
  }
  if (any(main_unassigned)) {
    cl_label <- match_n_assign(cl_label, which(assigned), main_unassigned, main_radius, distances)
  }
  cl_label
}


assign_unassigned <- function(distances,
                              cl_label,
                              seeds,
                              nng,
                              main_unassigned_method,
                              main_radius,
                              main_data_points,
                              secondary_unassigned_method,
                              secondary_radius) {
  assigned <- !is.na(cl_label)
  main_unassigned <- !assigned
  second_unassigned <- NULL
  if (!is.null(main_data_points)) {
    second_unassigned <- main_unassigned & !main_data_points
    main_unassigned <- main_unassigned & main_data_points
  }

  if (any(main_unassigned)) {
    if (main_unassigned_method == "ignore") {
      # nothing
    } else if (main_unassigned_method == "by_nng" || main_unassigned_method == "closest_assigned") {
      cl_label <- by_nng_or_closest_assigned(main_unassigned_method,
                                             main_unassigned,
                                             nng,
                                             assigned,
                                             main_radius,
                                             cl_label,
                                             distances)
    } else if (main_unassigned_method == "closest_seed") {
      cl_label <- match_n_assign(cl_label, seeds, main_unassigned, main_radius, distances)
    } else {
      stop("Unknown options.")
    }
  }

  if (!is.null(second_unassigned) && any(second_unassigned)) {
    if (secondary_unassigned_method == "ignore") {
      # nothing
    } else if (secondary_unassigned_method == "closest_assigned") {
      cl_label <- match_n_assign(cl_label, which(assigned), second_unassigned, secondary_radius, distances)
    } else if (secondary_unassigned_method == "closest_seed") {
      cl_label <- match_n_assign(cl_label, seeds, second_unassigned, secondary_radius, distances)
    } else {
      stop("Unknown options.")
    }
  }

  cl_label
}


est_average_seed_dist <- function(distances,
                                  nng,
                                  seeds) {
  if (length(seeds) > 1000) stop("Too many seeds for replica.")
  diag(nng) <- FALSE
  mean(unlist(lapply(seeds, function(i) {
    mean(distances[i, nng[, i]])
  })))
}


replica_nng_clustering <- function(distance_object,
                                   size_constraint,
                                   seed_method = "exclusion_updating",
                                   main_unassigned_method = "closest_seed",
                                   main_radius = NULL,
                                   main_data_points = NULL,
                                   secondary_unassigned_method = "ignore",
                                   secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  main_unassigned_method <- coerce_args(main_unassigned_method,
                                        c("ignore",
                                          "by_nng",
                                          "closest_assigned",
                                          "closest_seed",
                                          "estimated_radius_closest_seed"))
  main_radius <- coerce_radius(main_radius)
  if (is.null(main_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(main_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  distances <- as.matrix(distance_object)
  nng <- get_simple_nng(distances,
                        size_constraint,
                        main_radius,
                        main_data_points)
  seeds <- findseeds(nng, seed_method)
  cl_label <- assign_by_seed(seeds, nng)

  if (main_unassigned_method == "estimated_radius_closest_seed") {
    main_unassigned_method <- "closest_seed"
    main_radius <- est_average_seed_dist(distances, nng, seeds)
  }
  if (secondary_unassigned_method == "estimated_radius_closest_seed") {
    secondary_unassigned_method <- "closest_seed"
    secondary_radius <- est_average_seed_dist(distances, nng, seeds)
  }

  cl_label <- assign_unassigned(distances,
                                cl_label,
                                seeds,
                                nng,
                                main_unassigned_method,
                                main_radius,
                                main_data_points,
                                secondary_unassigned_method,
                                secondary_radius)

  make_Rscc_clustering(as.integer(cl_label),
                       length(unique(cl_label[!is.na(cl_label)])),
                       attr(distance_object, "ids", exact = TRUE))
}



replica_nng_clustering_types <- function(distance_object,
                                         type_labels,
                                         type_size_constraints,
                                         total_size_constraint = NULL,
                                         seed_method = "exclusion_updating",
                                         main_unassigned_method = "closest_seed",
                                         main_radius = NULL,
                                         main_data_points = NULL,
                                         secondary_unassigned_method = "ignore",
                                         secondary_radius = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count_distances(distance_object)
  type_labels <- coerce_type_labels(type_labels, num_data_points)
  type_size_constraints <- coerce_type_constraints(type_size_constraints)
  #ensure_type_labels_exist(names(type_size_constraints), get_all_types(type_labels))
  type_size_constraints <- make_type_size_constraints(type_size_constraints,
                                                      type_labels)
  total_size_constraint <- coerce_total_size_constraint(total_size_constraint,
                                                        type_size_constraints,
                                                        num_data_points)
  seed_method <- coerce_args(seed_method, all_seed_methods)
  main_unassigned_method <- coerce_args(main_unassigned_method,
                                        c("ignore",
                                          "by_nng",
                                          "closest_assigned",
                                          "closest_seed",
                                          "estimated_radius_closest_seed"))
  main_radius <- coerce_radius(main_radius)
  if (is.null(main_data_points)) {
    secondary_unassigned_method <- "ignore"
  } else {
    ensure_indicators(main_data_points, num_data_points, TRUE)
  }
  secondary_unassigned_method <- coerce_args(secondary_unassigned_method,
                                             c("ignore",
                                               "closest_assigned",
                                               "closest_seed",
                                               "estimated_radius_closest_seed"))
  secondary_radius <- coerce_radius(secondary_radius)

  distances <- as.matrix(distance_object)
  nng <- get_type_nng(distances,
                      type_labels,
                      type_size_constraints,
                      total_size_constraint,
                      main_radius,
                      main_data_points)
  seeds <- findseeds(nng, seed_method)
  cl_label <- assign_by_seed(seeds, nng)

  if (main_unassigned_method == "estimated_radius_closest_seed") {
    main_unassigned_method <- "closest_seed"
    main_radius <- est_average_seed_dist(distances, nng, seeds)
  }
  if (secondary_unassigned_method == "estimated_radius_closest_seed") {
    secondary_unassigned_method <- "closest_seed"
    secondary_radius <- est_average_seed_dist(distances, nng, seeds)
  }

  cl_label <- assign_unassigned(distances,
                                cl_label,
                                seeds,
                                nng,
                                main_unassigned_method,
                                main_radius,
                                main_data_points,
                                secondary_unassigned_method,
                                secondary_radius)

  make_Rscc_clustering(as.integer(cl_label),
                       length(unique(cl_label[!is.na(cl_label)])),
                       attr(distance_object, "ids", exact = TRUE))
}
