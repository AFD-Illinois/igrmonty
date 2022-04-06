project (igrmonty)

cmake_minimum_required(VERSION 3.19)
enable_testing()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
# Flags needed for GSL
set (LINK_FLAGS "-L/usr/local/lib -lgsl -lgslcblas -lm")
set (CMAKE_EXE_LINKER_FLAGS ${LINK_FLAGS})

configure_file(src/gff_ei_bremss.dat gff_ei_bremss.dat COPYONLY)
configure_file(src/bremss.c bremss.c COPYONLY)
configure_file(test/bremss/stub-decs.h decs.h COPYONLY)

add_executable(bremss_ei test/bremss/electron-ion.c)
add_test(NAME bremss_electron_ion COMMAND bremss_ei)

add_executable(bremss_ee test/bremss/electron-electron.c)
add_test(NAME bremss_electron_electron COMMAND bremss_ee)
