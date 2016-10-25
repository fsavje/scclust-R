if (identical(Sys.getenv("TRAVIS"), "true")) {
  run_slow_tests <- TRUE
} else if (identical(Sys.getenv("APPVEYOR"), "True")) {
  run_slow_tests <- TRUE
} else {
  run_slow_tests <- FALSE
}

compiled_with_stable_findseed <- FALSE

compiled_with_stable_nng <- FALSE
