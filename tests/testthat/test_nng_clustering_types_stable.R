library(Rscclust)
context("nng_clustering.R (types)")

source("../replica/replica_findseed.R", local = TRUE)
source("../replica/replica_make_nng.R", local = TRUE)
source("../replica/replica_nng.R", local = TRUE)
source("utils_nng.R", local = TRUE)


by_nng_or_closest_assigned <- function(main_unassigned_method,
                                       main_unassigned,
                                       nng,
                                       assigned,
                                       main_radius,
                                       cl_label,
                                       distances) {
  if (main_unassigned_method == "by_nng") {
    for(i in which(main_unassigned)) {
      pick <- intersect(which(nng[, i]), which(assigned))
      if (length(pick) > 0) {
        cl_label[i] <- cl_label[pick[1]]
        main_unassigned[i] <- FALSE
      }
    }
  } else if (any(main_unassigned)) {
    cl_label <- match_n_assign(cl_label, which(assigned), main_unassigned, main_radius, distances)
  }

  cl_label
}


test_that("type nng clustering function cluster correctly all combinations", {
  skip("Only run this when scclust is compiled with the -DSCC_STABLE_NNG flag.")

  for (tc in list(list(tc = "tc1",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 1L, "2" = 1L, "3" = 1L),
                       total_size_constraint = NULL),
                  list(tc = "tc2",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 2L, "2" = 2L, "3" = 1L),
                       total_size_constraint = 8L),
                  list(tc = "tc3",
                       type_labels = types1,
                       type_size_constraints = c("0" = 3L, "1" = 3L, "2" = 3L, "3" = 3L),
                       total_size_constraint = 12L),
                  list(tc = "tc4",
                       type_labels = types1,
                       type_size_constraints = c("0" = 1L, "1" = 0L, "2" = 1L, "3" = 0L),
                       total_size_constraint = 4L),
                  list(tc = "tc5",
                       type_labels = types2,
                       type_size_constraints = c("0" = 1L, "1" = 1L),
                       total_size_constraint = 2L),
                  list(tc = "tc6",
                       type_labels = types2,
                       type_size_constraints = c("0" = 2L, "1" = 0L),
                       total_size_constraint = 4L),
                  list(tc = "tc7",
                       type_labels = types2,
                       type_size_constraints = c("0" = 2L, "1" = 2L),
                       total_size_constraint = NULL))) {

    for (seed_method in c("lexical",
                          "inwards_order",
                          "inwards_updating",
                          "inwards_alt_updating",
                          "exclusion_order",
                          "exclusion_updating")) {

      for (main_unassigned_method in c("ignore",
                                       "by_nng",
                                       "closest_assigned",
                                       "closest_seed",
                                       "estimated_radius_closest_seed")) {

        radius_to_use <- c(0, 0.3, 0.4)
        if (isTRUE(tc$total_size_constraint > 6)) {
          radius_to_use <- c(0, 0.8)
        }
        if (isTRUE(tc$total_size_constraint > 10)) {
          radius_to_use <- c(0, 1.2)
        }

        for (main_radius in radius_to_use) {
          for (secondary_unassigned_method in c("ignore",
                                                "ignore_withmain",
                                                "closest_assigned",
                                                "closest_seed",
                                                "estimated_radius_closest_seed")) {
            for (secondary_radius in radius_to_use) {

              use_main_data_points <- NULL
              if (secondary_unassigned_method != "ignore") use_main_data_points <- main_data_points

              use_secondary_unassigned_method <- secondary_unassigned_method
              if (use_secondary_unassigned_method == "ignore_withmain") use_secondary_unassigned_method <- "ignore"

              use_main_radius <- main_radius
              if (use_main_radius == 0) use_main_radius <- NULL
              use_secondary_radius <- secondary_radius
              if (use_secondary_radius == 0) use_secondary_radius <- NULL

              cat(c(tc$tc,
                    seed_method,
                    main_unassigned_method,
                    use_main_radius,
                    is.null(use_main_data_points),
                    use_secondary_unassigned_method,
                    use_secondary_radius), "\n")
              test_nng_types_against_replica(test_distances1,
                                             tc$type_labels,
                                             tc$type_size_constraints,
                                             tc$total_size_constraint,
                                             seed_method,
                                             main_unassigned_method,
                                             use_main_radius,
                                             use_main_data_points,
                                             use_secondary_unassigned_method,
                                             use_secondary_radius)
            }
          }
        }
      }
    }
  }
})
