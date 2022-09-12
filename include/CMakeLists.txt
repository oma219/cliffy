set(MS_SOURCES  ms_rle_string.hpp doc_queries.hpp)
add_library(ms OBJECT ${MS_SOURCES})
target_link_libraries(ms common sdsl)
target_include_directories(ms PUBLIC "${r-index_SOURCE_DIR}/internal")

set(COMMON_SOURCES common.hpp)
add_library(common OBJECT ${COMMON_SOURCES})
target_link_libraries(common)

set(PFP_SOURCES dictionary.hpp parse.hpp pfp.hpp)
add_library(pfp OBJECT ${PFP_SOURCES})
target_link_libraries(pfp common sdsl divsufsort divsufsort64)