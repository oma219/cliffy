# add_executable(pfp_doc pfp_doc.cpp ref_builder.cpp)
# target_link_libraries(pfp_doc common pfp gsacak sdsl ri "-fopenmp")
# target_include_directories(pfp_doc PUBLIC "../include/")
# target_compile_options(pfp_doc PUBLIC "-std=c++17" "-fopenmp" "-mavx2" "-mavx512bw" "-mavx512f" 
# "-mavx512cd" "-mavx512dq" "-mavx512vl" "-mavx512ifma" "-mavx512vbmi" "-mavx512pf" )

add_executable(cliffy pfp_doc.cpp ref_builder.cpp minimizer_digest.cpp)
target_link_libraries(cliffy common pfp gsacak64 sdsl ri "-fopenmp")
target_include_directories(cliffy PUBLIC "../include/")
target_compile_options(cliffy PUBLIC "-std=c++17" "-fopenmp" "-DM64" "-march=native")

add_executable(mtest minimizer_digest_test.cpp minimizer_digest.cpp)
target_include_directories(mtest PUBLIC "../include")
target_compile_options(mtest PUBLIC "-std=c++17")



