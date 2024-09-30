## Submission

This is an update to fix the following issues:

* Removing specification of CXX11
* Updating naming in DESCRIPTION
* Replacing "class(x) == ..." with a better design pattern.

## Test environments 
* local Fedora install, R 4.4.1
* ubuntu (latest), R release, old-rel1 (on github)
* Microsoft Windows latest, R release (on github)
* Mac OS (latest), R release (on github)
* win-builder.r-project.org, R devel


## R CMD check results

  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wp,-D_GLIBCXX_ASSERTIONS’
    ‘-Wp,-U_FORTIFY_SOURCE,-D_FORTIFY_SOURCE=3’ ‘-march=x86-64’
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 1 note ✖

## Notes

This is a minor update that does not change functionality.
The note above is a local configuration issue and does not
show up on other builds.
