cmake_minimum_required(VERSION 3.8)
project(
    iterative_learning_mpc
    VERSION 0.1.0
    DESCRIPTION "Iterative-Learning Model Predictive Control"
    LANGUAGES CXX
)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()


# Set default C11
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
  set(CMAKE_C_STANDARD_REQUIRED True)
endif()

# Set default C++23
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 23)
  set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()


################
# DEPENDENCIES #
################

find_package(Eigen3 REQUIRED)

find_package(GTest REQUIRED)

find_package(utils REQUIRED)


#############
# LIBRARIES #
#############

# # iterative_learning_mpc
# add_library(iterative_learning_mpc STATIC src/iterative_learning_mpc.cpp)
# target_include_directories(iterative_learning_mpc PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
# )
# target_link_libraries(iterative_learning_mpc PUBLIC Eigen3::Eigen utils::utils)
# set_target_properties(iterative_learning_mpc PROPERTIES POSITION_INDEPENDENT_CODE ON)


###############
# EXECUTABLES #
###############

# main
# add_executable(main src/main.cpp)
# target_include_directories(src/xmain PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
# )
# target_link_libraries(main PRIVATE iterative_learning_mpc)


#########
# TESTS #
#########

# # Enable testing with CTest
# enable_testing()

# # test_iterative_learning_mpc
# add_executable(test_iterative_learning_mpc tests/test_iterative_learning_mpc.cpp)
# target_include_directories(test_iterative_learning_mpc PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
# )
# target_link_libraries(test_iterative_learning_mpc PRIVATE
#   iterative_learning_mpc
#   GTest::gmock_main
#   GTest::gtest
# )


# # Include gtest
# include(GoogleTest)
# gtest_discover_tests(test_iterative_learning_mpc)





##############################################################################
#                                    ROS2                                    #
##############################################################################

if(DEFINED ENV{ROS_VERSION})

  ################
  # DEPENDENCIES #
  ################
  find_package(ament_cmake REQUIRED)

  find_package(rclcpp REQUIRED)
  find_package(tf2 REQUIRED)
  find_package(tf2_ros REQUIRED)

  find_package(nav_msgs REQUIRED)
  find_package(visualization_msgs REQUIRED)
  find_package(tf2_geometry_msgs REQUIRED)


  #############
  # C++ NODES #
  #############

  # iterative_learning_mpc_node
  # add_executable(iterative_learning_mpc_node src/iterative_learning_mpc_node.cpp)
  # target_include_directories(iterative_learning_mpc_node PUBLIC
  #   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  #   $<INSTALL_INTERFACE:include>
  # )
  # target_link_libraries(iterative_learning_mpc_node iterative_learning_mpc)
  # ament_target_dependencies(iterative_learning_mpc_node
  #   rclcpp
  #   tf2
  #   tf2_ros
  #   nav_msgs
  #   visualization_msgs
  #   tf2_geometry_msgs
  # )
  # target_compile_features(iterative_learning_mpc_node PUBLIC c_std_11 cxx_std_23)


  ##############
  # INSTALLING #
  ##############

  # # installing executables
  # install(
  #   TARGETS
  #     iterative_learning_mpc_node
  #   DESTINATION
  #   lib/${PROJECT_NAME}
	# )

  # # Installing directories
  # install(
	# 	DIRECTORY
	# 		config
	# 		launch
	# 	DESTINATION
	# 		share/${PROJECT_NAME}
  # )


  ###########
  # TESTING #
  ###########
  if(BUILD_TESTING)
    find_package(ament_lint_auto REQUIRED)
    # the following line skips the linter which checks for copyrights
    # comment the line when a copyright and license is added to all source files
    set(ament_cmake_copyright_FOUND TRUE)
    # the following line skips cpplint (only works in a git repo)
    # comment the line when this package is in a git repo and when
    # a copyright and license is added to all source files
    set(ament_cmake_cpplint_FOUND TRUE)
    ament_lint_auto_find_test_dependencies()
  endif()

  ament_package()
endif()