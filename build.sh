#!/usr/bin/bash

g++-7 -I $HOME/git_src/rapidjson/include -I commons check_virial_chi_kasusformat.cpp commons/*.cpp -O3 -g -o cstress
