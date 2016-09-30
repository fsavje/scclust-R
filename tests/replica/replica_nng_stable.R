by_nng_or_closest_assigned <- function(main_unassigned_method,
                                       main_unassigned,
                                       nng,
                                       assigned,
                                       main_radius,
                                       cl_label,
                                       distances) {
  for(i in which(main_unassigned)) {
    pick <- intersect(which(nng[, i]), which(assigned))
    if (length(pick) > 0) {
      cl_label[i] <- cl_label[pick[1]]
      main_unassigned[i] <- FALSE
    }
  }
  if (any(main_unassigned) && main_unassigned_method == "closest_assigned") {
    cl_label <- match_n_assign(cl_label, which(assigned), main_unassigned, main_radius, distances)
  }
  cl_label
}
