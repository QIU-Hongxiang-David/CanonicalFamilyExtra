## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

I have tested this package on Windows (locally, with`devtools::check_win_devel()`, and with `rhub`), MacOS (with https://mac.r-project.org/macbuilder/submit.html) and Linux (with `rhub`). Some notes I have found:

- "unable to verify current time": This is likely a local issue with my machine.
- "Examples with CPU (user + system) or elapsed time > 10s": The examples can run except that they take a bit long.
- "checking CRAN incoming feasibility" "Maintainer: 'Hongxiang Qiu <david940408@gmail.com>'": This package is submitted by the maintainer.
