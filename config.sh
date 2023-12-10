#!/bin/bash

##### CREATE SYMBOLIC LINKS ########
# Examples directory
cd examples
ln -s ../src/soqcs/libpysoqcs.so
ln -s ../src/soqcs/__init__.py  soqcs.py
cd ..

# Documentation links
ln -s ./doc/cpp_doc/html/index.html ./SOQCS_cpp.html
ln -s ./doc/phy_doc/build/html/index.html ./SOQCS_phy.html


