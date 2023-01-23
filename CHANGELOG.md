# PFP-Doc Change Log

## v1.0.2 - latest
- Updated the construction code to write out lcp values in 16-bit registers opposed to using 8 bit registers. This
  now allows values to be up to 2^16-1, so any values larger than that are rounded down to 2^16-1. 
- Updated the query subcommand to load the profiles using the 16-bit registers.
- Updated the SIMD code used for the predecessor table to use 16-bit integers, I had to use
  masked versions of various functions since they were not present in header.

## v1.0.1
- Added info subcommand for printing out a section of the document array profiles
- Updated gsacak repo to avoid error with compilation

## v1.0.0
- Initial version of the document array profile construction algorithm.
- Initial version of the querying method for the document array profiles, where it takes reads and lists all the documents for each exact match it encounters while doing backward search.