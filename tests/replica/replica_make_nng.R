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

