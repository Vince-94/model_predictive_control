#!/bin/bash

./build/linear_mpc/test_linear_mpc --gtest_color=yes --gtest_print_time=1
./build/linear_time_varying_mpc/test_linear_time_varying_mpc --gtest_color=yes --gtest_print_time=1
