## Resubmission

This is a resubmission. In this version I have:

 * revised the Author fields to note other sources of code and 
   copyright holders to that code.
 * fixed documentation of internal, non-exported functions so that
   no un-runnable examples exist.
 * added arXiv link to the methods implemented in the package.

This should address the comments on the previous submission.

## Test environments (this resubmission)
* local Fedora 30, R 3.6.0
* ubuntu 14.04 (on travis-ci), rel, devel, oldrel
* winbuilder (devel)


## R CMD check results

0 errors | 0 warnings | 1 notes

Some of the installations had notes about this 

* This is a new release.
* Found the following (possibly) invalid URLs:
  URL: http://www.jstor.org/stable/2337067
    From: inst/doc/matrix-t-estimation.html
    Status: 403
		Message: Forbidden
		
