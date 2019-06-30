## Test environments
* local Fedora 30, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0, devel, oldrel
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

The various installations had, between them, a few notes (inconsistently)

* This is a new release.
* Found the following (possibly) invalid URLs:
  URL: http://www.jstor.org/stable/2337067
    From: inst/doc/matrix-t-estimation.html
    Status: 403
    Message: Forbidden
	
	This URL is generated from the valid DOI and works.
*  Found the following \keyword or \concept entries
   which likely give several index terms:
   File ‘rmatrixt.Rd’:
    \concept{matrix variate distributions

	I'm using \concept{matrix variate distributions} because it seems to be an 
	accurate description, I can use a different one if that causes issues.
