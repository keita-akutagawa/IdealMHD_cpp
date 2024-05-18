#!/bin/bash

options="-lgtest -lgtest_main -pthread"
file="test_calculate_half_components.cpp ../muscl.cpp ../minmod.cpp"
g++ $file $options
./a.out


