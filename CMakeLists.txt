cmake_minimum_required(VERSION 3.20)

project(
    model_predictive_control
    VERSION 0.1.0
    DESCRIPTION "Model Predictive Control stack"
    LANGUAGES CXX
)


##################
# GLOBAL OPTIONS #
##################

# Set default C++
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 23)
  set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

# Define compiler options
set(GCC_CLANG_FLAGS "-Wall" "-Wextra" "-Wpedantic" "-Werror")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "${CMAKE_CXX_COMPILER_ID} compiler detected")
  add_compile_options(${GCC_CLANG_FLAGS})
else()
  message(STATUS "Compiler is ${CMAKE_CXX_COMPILER_ID}")
endif()

# Set default build type
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

# Build type (default to Release)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Warnings (per-target, not global)
set(PROJECT_WARNINGS
  -Wall
  -Wextra
  -Wpedantic
)

################
# DEPENDENCIES #
################
find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(GTest REQUIRED)

# find_package(proxsuite REQUIRED)
# find_package(yaml-cpp REQUIRED)


###########
# MODULES #
###########

include_directories(${PROJECT_SOURCE_DIR}/model_predictive_control/src/utils/include)

add_subdirectory(model_predictive_control/src/linear_mpc)
add_subdirectory(model_predictive_control/src/linear_time_varying_mpc)


#########
# TESTS #
#########
# enable_testing()

# add_executable(test_data_structure model_prediction_control/utils/tests/test_data_structure.cpp)
# target_include_directories(test_data_structure PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
# )
# target_link_libraries(test_data_structure PRIVATE
#   linear_mpc
#   GTest::gmock_main
#   GTest::gtest
# )

# include(GoogleTest)
# gtest_discover_tests(test_data_structure)
