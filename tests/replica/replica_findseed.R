findseeds <- function(nng,
                      seed_method) {
  seeds <- NULL
  if (seed_method == "lexical") {
    seeds <- findseeds_lexical(nng)
  } else if (seed_method == "inwards_order") {
    seeds <- findseeds_inwards_order(nng)
  } else if (seed_method == "inwards_updating") {
    seeds <- findseeds_inwards_updating(nng)
  } else if (seed_method == "inwards_alt_updating") {
    seeds <- findseeds_inwards_alt_updating(nng)
  } else if (seed_method == "exclusion_order") {
    seeds <- findseeds_exclusion_order(nng)
  } else if (seed_method == "exclusion_updating") {
    seeds <- findseeds_exclusion_updating(nng)
  }
  seeds
}


findseeds_byorder <- function(orderind,
                              nng) {
  seeds <- NULL
  assigned <- rep(FALSE, ncol(nng))
  invisible(lapply(orderind, function(i) {
    if (any(nng[, i]) && !any(assigned[nng[, i]])) {
      seeds <<- c(seeds, i)
      assigned[nng[, i]] <<- TRUE
    }
  }))

  seeds
}


findseeds_lexical <- function(nng) {
  findseeds_byorder(1:ncol(nng), nng)
}


findseeds_inwards_order <- function(nng) {
  findseeds_byorder(order(rowSums(nng)), nng)
}


findseeds_inwards_updating <- function(nng) {
  seeds <- NULL
  assigned <- rep(FALSE, ncol(nng))
  inwards_count <- rowSums(nng - diag(diag(nng)))
  check_order <- order(inwards_count)

  while (length(check_order) > 0) {
    i <- check_order[1]
    check_order <- check_order[-1]
    if (any(nng[, i]) && !any(assigned[nng[, i]])) {
      seeds <- c(seeds, i)
      assigned[nng[, i]] <- TRUE

      for (e in which(nng[, i])) {
        for (nne in intersect(which(nng[, e]), check_order)) {
          if (any(nng[, nne]) && !assigned[nne]) {
            swap_pos <- which(inwards_count[check_order] == inwards_count[nne])[1]
            check_order[check_order == nne] <- check_order[swap_pos]
            check_order[swap_pos] <- nne
            inwards_count[nne] <- inwards_count[nne] - 1
          }
        }
      }
    }
  }
  seeds
}


findseeds_inwards_alt_updating <- function(nng) {
  seeds <- NULL
  assigned <- rep(FALSE, ncol(nng))
  inwards_count <- rowSums(nng - diag(diag(nng)))
  check_order <- order(inwards_count)

  while (length(check_order) > 0) {
    i <- check_order[1]
    check_order <- check_order[-1]
    if (any(nng[, i]) && !any(assigned[nng[, i]])) {
      seeds <- c(seeds, i)
      assigned[nng[, i]] <- TRUE
      for (e in intersect(which(nng[, i]), check_order)) {
        for (nne in intersect(which(nng[, e]), check_order)) {
          if (any(nng[, nne]) && !assigned[nne]) {
            swap_pos <- which(inwards_count[check_order] == inwards_count[nne])[1]
            check_order[check_order == nne] <- check_order[swap_pos]
            check_order[swap_pos] <- nne
            inwards_count[nne] <- inwards_count[nne] - 1
          }
        }
      }
    } else if (!assigned[i]) {
      for (nne in intersect(which(nng[, i]), check_order)) {
        if (any(nng[, nne]) && !assigned[nne]) {
          swap_pos <- which(inwards_count[check_order] == inwards_count[nne])[1]
          check_order[check_order == nne] <- check_order[swap_pos]
          check_order[swap_pos] <- nne
          inwards_count[nne] <- inwards_count[nne] - 1
        }
      }
    }
  }
  seeds
}


findseeds_exclusion_order <- function(nng) {
  exclusion_graph <- ((nng + t(nng) + t(nng) %*% nng) > 0)
  exclusion_graph[!apply(nng, 2, any), ] <- FALSE
  diag(exclusion_graph) <- FALSE
  findseeds_byorder(order(colSums(exclusion_graph)), nng)
}


findseeds_exclusion_updating <- function(nng) {
  seeds <- NULL
  excluded <- !apply(nng, 2, any)
  exclusion_list <- lapply(1:ncol(nng), function(i) {
    out <- c(which(nng[, i]), rev(which(nng[i, ])))
    for(e in which(nng[, i])) {
      out <- c(out, rev(which(nng[e, ])))
    }
    out <- out[!excluded[out]]
    out <- setdiff(out, i)
    unique(out)
  })
  exclud_count <- sapply(exclusion_list, length)
  check_order <- order(exclud_count)

  while (length(check_order) > 0) {
    i <- check_order[1]
    check_order <- check_order[-1]
    if (!excluded[i]) {
      seeds <- c(seeds, i)
      excluded[i] <- TRUE
      to_exclude <- intersect(exclusion_list[[i]], which(!excluded))
      excluded[to_exclude] <- TRUE
      for (e in to_exclude) {
        for (nne in intersect(exclusion_list[[e]], which(!excluded))) {
          swap_pos <- which(exclud_count[check_order] == exclud_count[nne])[1]
          check_order[check_order == nne] <- check_order[swap_pos]
          check_order[swap_pos] <- nne
          exclud_count[nne] <- exclud_count[nne] - 1
        }
      }
    }
  }
  seeds
}
