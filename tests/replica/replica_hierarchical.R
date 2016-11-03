find_centers <- function(indices, distances) {
  ISCC_GR_NUM_TO_CHECK <- 100L
  step <- max(2L, length(indices) %/% ISCC_GR_NUM_TO_CHECK)
  to_check <- indices[seq(1, min(step * ISCC_GR_NUM_TO_CHECK, length(indices)), step)]
  checked <- to_check
  best <- NULL

  while (length(to_check) > 0) {
    check_dist <- distances[to_check, indices, drop = FALSE]
    candidate <- arrayInd(which.max(t(check_dist)), dim(t(check_dist)))
    cand_tmp <- as.integer(to_check[candidate[1,2]])
    candidate[1,2] <- as.integer(indices[candidate[1,1]])
    candidate[1,1] <- cand_tmp
    if (is.null(best) || (distances[candidate] > distances[best])) {
      best <- candidate
    }
    to_check <- setdiff(indices[unique(apply(check_dist, 1, which.max))], checked)
    checked <- c(checked, to_check)
  }

  best
}


break_cluster <- function(indices,
                          distances,
                          size_constraint,
                          batch_assign) {
  centers <- as.integer(find_centers(indices, distances))
  clusters <- list(centers[1], centers[2])
  unassigned <- setdiff(indices, centers)

  if (sort(distances[centers[1], unassigned])[size_constraint - 1] >=
      sort(distances[centers[2], unassigned])[size_constraint - 1]) {
    clusters[[1]] <- c(clusters[[1]], unassigned[order(distances[centers[1], unassigned])[1:(size_constraint - 1)]])
    unassigned <- setdiff(unassigned, clusters[[1]])
    clusters[[2]] <- c(clusters[[2]], unassigned[order(distances[centers[2], unassigned])[1:(size_constraint - 1)]])
    unassigned <- setdiff(unassigned, clusters[[2]])
  } else {
    clusters[[2]] <- c(clusters[[2]], unassigned[order(distances[centers[2], unassigned])[1:(size_constraint - 1)]])
    unassigned <- setdiff(unassigned, clusters[[2]])
    clusters[[1]] <- c(clusters[[1]], unassigned[order(distances[centers[1], unassigned])[1:(size_constraint - 1)]])
    unassigned <- setdiff(unassigned, clusters[[1]])
  }

  while (length(unassigned) > 0) {
    if (batch_assign) {
      if (length(unassigned) < size_constraint) {
        size_constraint <- length(unassigned)
      }
      if (sort(distances[centers[1], unassigned])[size_constraint] <=
          sort(distances[centers[2], unassigned])[size_constraint]) {
        clusters[[1]] <- c(clusters[[1]], unassigned[order(distances[centers[1], unassigned])[1:size_constraint]])
        unassigned <- setdiff(unassigned, clusters[[1]])
      } else {
        clusters[[2]] <- c(clusters[[2]], unassigned[order(distances[centers[2], unassigned])[1:size_constraint]])
        unassigned <- setdiff(unassigned, clusters[[2]])
      }
    } else {
      if (sort(distances[centers[1], unassigned])[1] <=
          sort(distances[centers[2], unassigned])[1]) {
        clusters[[1]] <- c(clusters[[1]], unassigned[order(distances[centers[1], unassigned])[1]])
        unassigned <- setdiff(unassigned, clusters[[1]])
      } else {
        clusters[[2]] <- c(clusters[[2]], unassigned[order(distances[centers[2], unassigned])[1]])
        unassigned <- setdiff(unassigned, clusters[[2]])
      }
    }
  }

  clusters[[2]] <- rev(clusters[[2]])

  clusters
}


run_hierarchical <- function (num_data_points,
                              cluster_queue,
                              distances,
                              size_constraint,
                              batch_assign) {
  cl_labels <- rep(NA, num_data_points)
  current_label <- 0L
  cluster_queue <- rev(cluster_queue)

  while(length(cluster_queue) > 0) {
    if (length(cluster_queue[[1]]) < 2 * size_constraint) {
      cl_labels[cluster_queue[[1]]] <- current_label
      current_label <- current_label + 1L
      cluster_queue <- cluster_queue[-1]
    } else {
      new_clusters <- break_cluster(cluster_queue[[1]], distances, size_constraint, batch_assign)
      cluster_queue[[1]] <- new_clusters[[1]]
      cluster_queue <- c(new_clusters[2], cluster_queue)
    }
  }

  cl_labels
}


replica_hierarchical_clustering <- function(distance_object,
                                            size_constraint,
                                            batch_assign = TRUE,
                                            existing_clustering = NULL) {
  ensure_distances(distance_object)
  num_data_points <- data_point_count.Rscc_distances(distance_object)
  size_constraint <- coerce_size_constraint(size_constraint, num_data_points)
  ensure_indicators(batch_assign, 1L)
  if (!is.null(existing_clustering)) {
    ensure_Rscc_clustering(existing_clustering, num_data_points)
  }

  if (!is.null(existing_clustering)) {
    cluster_queue <- lapply(sort(unique(existing_clustering)), function(x) { rev(which(existing_clustering == x)) })
  } else {
    cluster_queue <- list(1:num_data_points)
  }

  new_labels <- run_hierarchical(num_data_points,
                                 cluster_queue,
                                 as.matrix(distance_object),
                                 size_constraint,
                                 batch_assign)

  make_Rscc_clustering(new_labels,
                       length(unique(new_labels)),
                       attr(distance_object, "ids", exact = TRUE))
}
