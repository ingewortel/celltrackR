## Resubmission
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


## Test environments
* local OS X install, R 3.6.1
* ubuntu 18.04.3, R 3.6.2
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Johannes Textor <johannes.textor@gmx.de>’

New submission
