findseeds_inwards_updating <- function(nng) {
  seeds <- NULL
  assigned <- rep(FALSE, ncol(nng))
  tobechecked <- 1:ncol(nng)

  while (length(tobechecked) > 0) {
    i <- tobechecked[order(rowSums(nng[tobechecked, intersect(which(!assigned), tobechecked), drop = FALSE]))[1]]
    tobechecked <- setdiff(tobechecked, i)
    if (any(nng[, i]) && !any(assigned[nng[, i]])) {
      seeds <- c(seeds, i)
      tobechecked <- setdiff(tobechecked, which(nng[, i]))
      assigned[nng[, i]] <- TRUE
    }
  }

  seeds
}


findseeds_exclusion_updating <- function(nng) {
  seeds <- NULL
  exclusion_graph <- ((nng + t(nng) + t(nng) %*% nng) > 0)
  exclusion_graph[!apply(nng, 2, any), ] <- FALSE
  exclud_count <- colSums(exclusion_graph)
  tobechecked <- which(apply(nng, 2, any))

  while (length(tobechecked) > 0) {
    i <- tobechecked[order(exclud_count[tobechecked])[1]]
    tobechecked <- setdiff(tobechecked, i)
    seeds <- c(seeds, i)
    to_exclude <- intersect(which(exclusion_graph[, i]), tobechecked)
    tobechecked <- setdiff(tobechecked, to_exclude)

    for (e in to_exclude) {
      for (nne in intersect(which(exclusion_graph[, e]), tobechecked)) {
        exclud_count[nne] <- exclud_count[nne] - 1
      }
    }
  }

  seeds
}
