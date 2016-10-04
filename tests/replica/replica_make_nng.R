get_simple_nng <- function(distances,
                           num_nn,
                           radius,
                           main_data_points) {
  num_data_points <- ncol(distances)
  nng <- matrix(FALSE, ncol = num_data_points, nrow = num_data_points)

  if (is.null(main_data_points)) {
    lookup <- 1:num_data_points
  } else {
    lookup <- which(main_data_points)
  }

  invisible(lapply(lookup, function(i) {
    nn <- order(distances[i, ])[1:num_nn]
    if (is.null(radius) || distances[i, nn[num_nn]] <= radius) {
      nng[nn, i] <<- TRUE
    }
  }))

  nng
}


get_type_nng <- function(distances,
                         type_labels,
                         type_size_constraints,
                         total_size_constraint,
                         main_radius,
                         main_data_points) {
  num_data_points <- ncol(distances)
  nng <- matrix(FALSE, ncol = num_data_points, nrow = num_data_points)

  if (is.null(main_data_points)) {
    lookup <- 1:num_data_points
  } else {
    lookup <- which(main_data_points)
  }

  for (t in 1:length(type_size_constraints)) {
    if (type_size_constraints[t] > 0) {
      new_lookup <- lookup
      of_type <- which(type_labels == (t - 1))
      for (i in lookup) {
        type_dist <- as.numeric(distances[i, of_type])
        type_dist_order <- order(type_dist)[1:type_size_constraints[t]]
        if (!is.null(main_radius) && type_dist[type_dist_order[type_size_constraints[t]]] > main_radius) {
          nng[, i] <- FALSE
          new_lookup <- setdiff(new_lookup, i)
        } else {
          nng[of_type[type_dist_order], i] <- TRUE
        }
      }
      lookup <- new_lookup
    }
  }

  extra_con <- total_size_constraint - sum(type_size_constraints)
  if (extra_con > 0) {
    for (i in lookup) {
      unconnected <- which(!nng[, i])
      type_dist <- as.numeric(distances[i, unconnected])
      type_dist_order <- order(type_dist)[1:extra_con]
      if (!is.null(main_radius) && type_dist[type_dist_order[extra_con]] > main_radius) {
        nng[, i] <- FALSE
      } else {
        nng[unconnected[type_dist_order], i] <- TRUE
      }
    }
  }

  nng
}
