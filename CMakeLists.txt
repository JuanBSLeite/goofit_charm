cmake_minimum_required(VERSION 3.5)

goofit_add_package(CHARM)

goofit_add_executable(D2PPP D2PPP.cpp)
goofit_add_executable(D2PPP_MIPWA D2PPP_MIPWA.cpp)
#if(GOOFIT_DEVICE STREQUAL CUDA)
#    target_compile_options(D2PPP PUBLIC --expt-extended-lambda)
#endif()
