# PFP-Doc Change Log

## v1.0.6 - latest
- Added another method to trim the construction queue in scenarios such as large block of suffixes that 
  begin with Ns or any character that is not a usual DNA/RNA character.

## v1.0.5
- Used the second method for lcp queue trimming which is a heuristic, anytime we encounter an lcp value <= 5, we flush the queue
- Made a change to update_lcp_queue to incorporate knowledge of which method was used for trimming

## v1.0.4 
- Turned off second method for lcp queue trimming to avoid errors
- Fixed the non-heuristic method for lcp queue trimming by updating the ch_doc_counter table 
  during the loop
- Updated the code that updates profile with LF steps to avoid overflow

## v1.0.3
- Updated the doc_queries constructor to avoid creating a separate table of sequential entries
  to limit the RAM usage, especially when working with large datasets.

## v1.0.2
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