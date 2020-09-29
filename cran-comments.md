## Resubmission

This version fixes issues caused by an update of the `jrc` package to version 0.4.0

## Test environments
* local ubuntu 20.04 LST, R 4.0.2
* win-builder: R-devel
* ubuntu 16.04.6 LTS, R-devel, R 4.0.2, R 3.6.3 (travis ci)

## R CMD check results

There were no ERRORS, WARNINGS or NOTES

# Previous cran-comments

## Resubmission

This is a new version of the package. It adds new argument `metric` that provides an alternative
to default Euclidean distance.

## Test environments
* local ubuntu 18.04 LST, R 3.6.2
* win-builder: R-devel
* ubuntu 16.04.6 LTS, R-devel, R 3.6.2, R 3.5.3 (travis ci)

## R CMD check results

There were no ERRORS, WARNINGS or NOTES

## Resubmission

This is a new version of the package. The changes are caused by update `jrc` package.

## Test environments
* local ubuntu 18.04 LST, R 3.6.1
* win-builder: R-devel
* ubuntu 14.04.5 LTS, R-devel, R 3.5.2, R 3.4.4 (travis ci)

## R CMD check results

There were no ERRORS, WARNINGS or NOTES

## Resubmission

This is a new version of the package. Some performance improvent was made and a minor bug fixed. See News.md
for more details.

## Test environments
* local ubuntu 18.04 LST, R 3.6.1
* win-builder: R-devel
* ubuntu 14.04.5 LTS, R-devel, R 3.5.2, R 3.4.4 (travis ci)

## R CMD check results

There were no ERRORS, WARNINGS or NOTES

## Test environments
* local ubuntu 16.04 LST, R 3.5.1
* win-builder: R-devel
* ubuntu 14.04.5 LTS, R-devel, R 3.5.2, R 3.4.4 (travis ci)

## R CMD check results

There were no ERRORS, WARNINGS or NOTES

## Resubmission

* DESCRIPTION changed

* slw_on_selection is no longer an exported function. Instead a 
	variable to store the apps state is exported (.slw).

* an HTML filed created in the example is now saved to tempdir(),
	sleepwalk function doesn't save any files by default.

* We've just submitted the related paper, and therefore no reference can be provided at the
	moment.

>> Please provide executable examples in slw_on_selection.Rd and unwrap /donttest in slw_snapshot.Rd
Unfortunately, slw_snapshot requires an esteblished WebSocket connection with an opened web page, 
which is impossible during the automated testing. If there is no currently running Sleepwalk page,
the function will wait for 10 seconds and produce an error "Failed to get embedding data from the server".

