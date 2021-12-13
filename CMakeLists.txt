# ======================================================================
#  larsimrad main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsimrad/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larsimrad VERSION 09.04.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

# -Wno-unused-variable needed because of CLHEP
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-unused-variable
)

cet_report_compiler_flags()

find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package(nugen REQUIRED PUBLIC)
find_package(nurandom REQUIRED PUBLIC)
find_package(GSL REQUIRED)

find_ups_product( larcore )
find_ups_product( larcoreobj )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larsim )
find_ups_product( bxdecay0 )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# source
add_subdirectory(larsimrad)

# tests
#add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
