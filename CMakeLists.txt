project (igrmonty)

cmake_minimum_required(VERSION 3.19)
enable_testing()

set(CMAKE_C_STANDARD 99)

# prereqs
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)
# TODO can/should we avoid this req't for non-iharm models?
find_package(HDF5 COMPONENTS C)

# includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/model/template)
include_directories(${GSL_INCLUDE_DIRS})

# TESTS
configure_file(src/gff_ei_bremss.dat gff_ei_bremss.dat COPYONLY)

add_executable(bremss_ei test/bremss/electron-ion.c)
target_link_libraries(bremss_ei GSL::gsl GSL::gslcblas)
target_link_libraries(bremss_ei m)
add_test(NAME bremss_electron_ion COMMAND bremss_ei)

add_executable(bremss_ee test/bremss/electron-electron.c)
target_link_libraries(bremss_ee GSL::gsl GSL::gslcblas)
target_link_libraries(bremss_ee m)
add_test(NAME bremss_electron_electron COMMAND bremss_ee)

# GRMONTY
# Spedify a model with MODEL=
set(MODEL "iharm" CACHE STRING "Model to compile")
# Extra includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/model/${MODEL})
# Enumerate source files, compile
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src EXE_NAME_SRC)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/model/${MODEL} EXE_NAME_SRC)
add_executable(grmonty ${EXE_NAME_SRC})
target_link_libraries(grmonty m)
target_link_libraries(grmonty GSL::gsl GSL::gslcblas)
target_link_libraries(grmonty OpenMP::OpenMP_CXX)
target_link_libraries(grmonty ${HDF5_C_LIBRARIES})