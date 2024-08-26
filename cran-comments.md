## Current submission v1.2.1

### Resubmission after archiving

The package was archived due to a violation of the policy on handling internet resources,
which we were unable to fix in time before the deadline (2024-08-20). We have now
adapted the package in line with the policy and are resubmitting it to CRAN. 


### Test environments
* local OS X (aarch64-apple-darwin20 (64-bit) with R v4.4.1 (2024-08-23) 
* win-builder (devel) (2024-08-23)
* windows-latest (release) x86_64-w64-mingw32 (64-bit) with R v4.4.1 (2024-08-23)
* ubuntu-latest (devel) x86_64-pc-linux-gnu (64-bit) with R developmental version (2024-08-23)
* ubuntu-latest (release) x86_64-pc-linux-gnu (64-bit) with R v4.4.1 (2024-08-23)


### R CMD check results

Status: OK

0 errors | 0 warnings | 0 notes

### Reverse dependencies

#### revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages



## Previous submission v1.2.0


### Test environments
* local OS X (x86_64-apple-darwin17.0 (64-bit) with R v4.2.3 (2024-03-25)
* win-builder (devel) (2024-03-25)
* windows-latest (release) x86_64-w64-mingw32 (64-bit) with R v4.3.2 (2024-03-25)
* ubuntu-latest (devel) x86_64-apple-darwin20 (64-bit) with R v4.3.2 (2024-03-25)
* ubuntu-latest (release) x86_64-pc-linux-gnu (64-bit) with R v4.3.2 (2024-03-25)


### R CMD check results

Status: OK

0 errors | 0 warnings | 0 notes

### Reverse dependencies

#### revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


## Previous submission v1.1.0

### Resubmission

Edit 22-02-2022: Upon request, we have:

* inst/CITATION : added trailing slash to URL as needed
* vignettes/*.Rmd, inst/doc/*.html : updated moved URLs as needed

### Test environments
* local OS X (x86_64-apple-darwin20.3.0, 64-bit) with R v4.0.5 (2021-03-31)
* ubuntu 20.04, R v3.6.3
* win-builder (devel)

### R CMD check results

Status: OK

(no errors/warnings/notes).

### Reverse dependencies

Not applicable: no dependencies found with revdepcheck::revdep_check().


## Older submission

### Resubmission
Edit 30-03-2020 (2): Upon request, we have:

* DESCRIPTION : updated the date field.

Edit 30-03-2020: Upon request, we have made the following, additional changes:

* DESCRIPTION : Added references as requested;
* R/ : Fixed problem with resetting par() settings twice in one function;
* man/ : Added \value to .Rd file TrackMeasures.Rd.

This is a resubmission. As requested, we have made the following changes in this version:

* DESCRIPTION : Updated the Author section using the Authors@R field, declaring 
	Maintainer, Authors and Contributors with their appropriate roles;
* R/ : Used on.exit() to reset par() settings on function exist for all functions that 
	change these settings;
* man/ : Added \value to all .Rd files of exported functions where this was missing (In the
	corresponding code in the R/ directory, this corresponds to adding a @return tag
	in the comments used by Roxygen);
* vignettes/ : Reset par() at the end of every .Rmd vignette;
* vignettes/ : Renamed any variables in the vignettes named T.[something];


### Test environments
* local OS X install, R 3.6.1
* ubuntu 18.04.3, R 3.6.2
* win-builder (devel)

### R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Johannes Textor <johannes.textor@gmx.de>’

New submission
