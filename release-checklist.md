## Submit release

* Check and fix any errors at https://cran.r-project.org/web/checks/check_results_scclust.html
* Delete content of `NAMESPACE` and run `devtools::document()`
* Run `devtools::load_all(recompile = TRUE)` and `devtools::check()`
* Run `devtools::check_win_devel()`, `devtools::check_win_release()` and `devtools::check_win_oldrelease()`
* Run `revdepcheck::revdep_check()`, remove "revdep" folder when done
* Run comprehensive tests locally
	- Change `run_slow_tests <- TRUE` in `tests/testthat/config.R`
	- Run `devtools::load_all(recompile = TRUE)` and `devtools::test()` with the following settings
	- Remember to change the flag in the Makefile as well
		- stable_nng = FALSE, stable_findseed = FALSE (medium time)
		- stable_nng = TRUE, stable_findseed = FALSE (long time)
		- stable_nng = FALSE, stable_findseed = TRUE (short time)
		- stable_nng = TRUE, stable_findseed = TRUE (long time)
	- Change `tests/testthat/config.R` and Makefile to original state
* Update package information
	- Set new version number in `DESCRIPTION`
	- Set release date in `DESCRIPTION`
	- Change "scclust devel" to "scclust VERSION" in `NEWS.md`
	- Update `cran-comments.md` with correct information
	- Update travis and appveyor with current versions
* Commit and push to github so automatic tests run
* Run `devtools::load_all(recompile = TRUE)`, `devtools::test()` and `devtools::check()`
* Run `devtools::check_win_devel()`, `devtools::check_win_release()` and `devtools::check_win_oldrelease()`
* Run `revdepcheck::revdep_check()`, remove "revdep" folder when done
* Wait until all tests are done
* Submit to CRAN
    - Run `devtools::build()`
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
