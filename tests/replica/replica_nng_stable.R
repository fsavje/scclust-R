by_nng_or_closest_assigned <- function(unassigned_method,
                                       unassigned,
                                       nng,
                                       assigned,
                                       radius,
                                       cl_label,
                                       distances) {
  for(i in which(unassigned)) {
    pick <- intersect(which(nng[, i]), which(assigned))
    if (length(pick) > 0) {
      cl_label[i] <- cl_label[pick[1]]
      unassigned[i] <- FALSE
    }
  }
  if (any(unassigned) && unassigned_method == "closest_assigned") {
    cl_label <- match_n_assign(cl_label, which(assigned), unassigned, radius, distances)
  }
  cl_label
}
