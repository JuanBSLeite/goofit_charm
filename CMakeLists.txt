cmake_minimum_required(VERSION 3.8)

goofit_add_package(CHARM)

goofit_add_executable(D2PPP D2PPP.cpp)

#if(GOOFIT_DEVICE STREQUAL CUDA)
#    target_compile_options(D2PPP PUBLIC --expt-extended-lambda)
#endif()
