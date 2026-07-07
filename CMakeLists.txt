cmake_minimum_required(VERSION 3.20)
project(sedpp CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -g")

find_package(PkgConfig REQUIRED)
pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)
pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)
pkg_check_modules(NLOHMANN_JSON REQUIRED IMPORTED_TARGET nlohmann_json)
pkg_check_modules(CURL REQUIRED IMPORTED_TARGET libcurl)
pkg_check_modules(ZLIB REQUIRED IMPORTED_TARGET zlib)  # gzip cache of downloaded spectra

add_library(mpfit STATIC extern/mpfit/mpfit.c)
target_include_directories(mpfit PUBLIC extern/mpfit)

add_executable(sedfit
  src/main.cpp
  src/config.cpp
  src/photometry_table.cpp
  src/passbands.cpp
  src/synthmag.cpp
  src/extinction.cpp
  src/modelgrid.cpp
  src/fitsio_util.cpp
  src/fitfun.cpp
  src/fitter.cpp
  src/stellar.cpp
  src/texval.cpp
  src/texout.cpp
  src/fitsout.cpp
  src/results.cpp
  src/magflux.cpp
  src/spectrum.cpp
  src/util.cpp
  src/http.cpp
  src/votable.cpp
  src/catalog_registry.cpp
  src/reddening.cpp
  src/query.cpp
  src/spectra.cpp
  src/astrometry.cpp
  src/simbad.cpp
)
target_link_libraries(sedfit PRIVATE mpfit PkgConfig::CFITSIO PkgConfig::GSL PkgConfig::NLOHMANN_JSON PkgConfig::CURL PkgConfig::ZLIB)
