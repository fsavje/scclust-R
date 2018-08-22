## Submit release

* Check and fix any errors at https://cran.r-project.org/web/checks/check_results_scclust.html
* Delete content of `NAMESPACE` and run `document()`
* Run `load_all(recompile = TRUE)` and `check()`
* Run `build_win(version = "R-release")` and `build_win(version = "R-devel")`
* Run `revdep_check()`, remove `revdep` folder when done
* Run comprehensive tests locally
	- Change `run_slow_tests <- TRUE` in `tests/testthat/config.R`
	- Run `load_all(recompile = TRUE)` and `test()` with the following settings
	- Remember to change the flag in the Makefile as well
		- stable_nng = FALSE, stable_findseed = FALSE
		- stable_nng = TRUE, stable_findseed = FALSE
		- stable_nng = FALSE, stable_findseed = TRUE
		- stable_nng = TRUE, stable_findseed = TRUE
	- Change `tests/testthat/config.R` and Makefile to original state
* Update package information
	- Set new version number in `DESCRIPTION`
	- Set release date in `DESCRIPTION`
	- Change "scclust devel" to "scclust VERSION" in `NEWS.md`
	- Update `cran-comments.md` with correct information
	- Update travis and appveyor with current versions
* Commit and push to github so automatic tests run
* Run `load_all(recompile = TRUE)`, `test()` and `check()`
* Run `build_win(version = "R-release")` and `build_win(version = "R-devel")`
* Run `revdep_check()`, remove revdep folder when done
* Wait until all tests are done
* Submit to CRAN
	- Run `build()`
	- Upload to http://cran.r-project.org/submit.html
	- Add `cran-comments.md` as comment


## When accepted

* Add new release to Github
	- Add CRAN release as binary
	- Add relevant information from `NEWS.md`
* Update package information
	- Add .9000 to the version number in `DESCRIPTION`
	- Set date in `DESCRIPTION`
	- Add "scclust devel" to `NEWS.md`
	- Commit and push to github
