-## Submit comment

One downstream dependency (quickmatch) has a failed test with the new version
of scclust. One of the fixes to scclust in this version is that an erroneous error
message is no longer generated. The test that quickmatch fails is exactly that it
expect that error message from scclust. Note that the error message was produced
without any error present. Thus, the user experience is unaffected by the failed
test in quickmatch, so no existing code should be affected. quickmatch will be
updated to have proper tests as soon as the new version of scclust is on CRAN.


## Test environments

  * x86_64-apple-darwin15.6.0 (local machine)
     - 3.5.1

  * x86_64-w64-mingw32 (win-builder)
     - 3.5.1, devel

  * x86_64-pc-linux-gnu (travis-ci)
     - 3.4.0, 3.4.1, 3.4.2, 3.4.3, 3.4.4, 3.5.0,
       3.5.1, devel

  * x86_64-w64-mingw32/x64 (appveyor)
     - 3.4.0, 3.4.1, 3.4.2, 3.4.3, 3.4.4, 3.5.0,
       3.5.1, devel

  * i386-w64-mingw32/i386 (appveyor)
     - 3.4.4, 3.5.1, devel


## R CMD check results

  * No ERRORs

  * No WARNINGs.

  * No NOTEs.


## Downstream dependencies

  All downstream dependencies have been checked.

  * quickblock: 0 errors | 0 warnings | 0 notes

  * quickmatch: 1 errors | 0 warnings | 0 notes
