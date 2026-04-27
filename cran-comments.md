## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

I have used `rhub` to test this package on Windows (also with `devtools::check_win_devel()`), MacOS (also with https://mac.r-project.org/macbuilder/submit.html) and Linux. Some notes I have found:

- "unable to verify current time": This is likely an issue with my own machine.
- "Examples with CPU (user + system) or elapsed time > 10s": The examples can run except that they take a bit long.
- "checking CRAN incoming feasibility" "Maintainer: 'Hongxiang Qiu <david940408@gmail.com>'": This package is submitted by the maintainer.
