cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    MrLavaLoba
    VERSION 0.1.0
    DESCRIPTION "Stochastic model for lava flows"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    MrLavaLoba_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    MrLavaLoba_lib ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(MrLavaLoba_lib PUBLIC cxx_std_20)

find_package(fmt REQUIRED)
target_link_libraries(MrLavaLoba_lib PRIVATE fmt::fmt)

find_package(Eigen3 REQUIRED)
target_link_libraries(MrLavaLoba_lib PRIVATE Eigen3::Eigen)

find_package(Boost REQUIRED)
target_link_libraries(MrLavaLoba_lib PRIVATE Boost::program_options)


find_package(shapelib REQUIRED)
target_link_libraries(MrLavaLoba_lib PRIVATE shapelib::shp)

# ---- Declare executable ----

add_executable(MrLavaLoba_exe source/main.cpp)
add_executable(MrLavaLoba::exe ALIAS MrLavaLoba_exe)

set_property(TARGET MrLavaLoba_exe PROPERTY OUTPUT_NAME MrLavaLoba)

target_compile_features(MrLavaLoba_exe PRIVATE cxx_std_20)

target_link_libraries(MrLavaLoba_exe PRIVATE MrLavaLoba_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT MrLavaLoba_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of MrLavaLoba"
  )
endif()

include(cmake/dev-mode.cmake)
