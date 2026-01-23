cmake_minimum_required(VERSION 3.8)
project(utils)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()


################
# DEPENDENCIES #
################

find_package(Eigen3 REQUIRED)
find_package(GTest REQUIRED)


#############
# LIBRARIES #
#############

add_library(utils INTERFACE)
target_include_directories(utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(utils INTERFACE Eigen3::Eigen)

# Install headers (preserving utils/ subfolder)
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION include/${PROJECT_NAME}
)

# Install and export the INTERFACE target
install(
  TARGETS utils
  EXPORT export_${PROJECT_NAME}
)


#########
# TESTS #
#########

# # Enable testing with CTest
# enable_testing()

# # test_data_structure
# add_executable(test_data_structure tests/test_data_structure.cpp)
# target_include_directories(test_data_structure PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
# )
# target_link_libraries(test_data_structure PRIVATE
#   data_structure
#   Eigen3::Eigen
#   GTest::gmock_main
#   GTest::gtest
# )


# # Include gtest
# include(GoogleTest)
# gtest_discover_tests(test_data_structure)



##############################################################################
#                                    ROS2                                    #
##############################################################################

if(DEFINED ENV{ROS_VERSION})

  # find dependencies
  find_package(ament_cmake REQUIRED)

  # Tell ament this is the exported target (modern & recommended way)
  ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)
  ament_export_include_directories(include/${PROJECT_NAME})


  ament_package()
endif()