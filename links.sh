#!/bin/bash

##### CREATE SYMBOLIC LINKS ########
# Python directory
cd py_src
ln -s ../src
cd ..

# Development directory
cd devel
ln -s ../py_src/libSOQCS.so
ln -s ../py_src/pysoqcs.py
ln -s ../py_src/pysoqcs.conf
cd ..

# Examples directory
cd examples
ln -s ../py_src/libSOQCS.so
ln -s ../py_src/pysoqcs.py
ln -s ../py_src/pysoqcs.conf
cd ..

# Documentation links
ln -s ./doc/cpp_doc/html/index.html ./SOQCS_cpp.html
ln -s ./doc/phy_doc/build/html/index.html ./SOQCS_phy.html




