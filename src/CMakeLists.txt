add_executable(pfp_doc pfp_doc.cpp ref_builder.cpp)
target_link_libraries(pfp_doc common pfp gsacak sdsl ri)
target_include_directories(pfp_doc PUBLIC "../include/")
target_compile_options(pfp_doc PUBLIC "-std=c++17")