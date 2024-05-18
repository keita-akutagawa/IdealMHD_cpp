#!/bin/bash

options="-lgtest -lgtest_main -pthread"
programfile="../calculate_half_components.cpp ../muscl.cpp ../minmod.cpp ../const.cpp"
testfile="test_calculate_half_components.cpp"
g++ $testfile $programfile $options
./a.out


