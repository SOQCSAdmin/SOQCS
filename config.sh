#!/bin/bash

############################
##### CONFIGURATION ########
############################
VER=3.4.0            #VERSION
DIR=./cache          #CACHE DIRECTORY
INSTALL=./src/Eigen  #EIGEN INSTALL DIRECTION (DON'T TOUCH)


#####################
##### SCRIPT ########
#####################
##### UNPACK DOCUMENTATION ########
unzip doc.zip

##### INSTALL EIGEN 3 ########
mkdir $DIR
wget -O $DIR/eigen.tar.gz  https://gitlab.com/libeigen/eigen/-/archive/$VER/eigen-$VER.tar.gz
tar xvf $DIR/eigen.tar.gz -C $DIR
cp -R $DIR/eigen-$VER/Eigen/* $INSTALL/

##### CREATE SYMBOLIC LINKS ########
cd py_src
ln -s ../src
cd ..
cd devel
ln -s ../py_src/libSOQCS.so
ln -s ../py_src/pysoqcs.py
cd ..
cd examples
ln -s ../py_src/libSOQCS.so
ln -s ../py_src/pysoqcs.py
cd ..
ln -s ./doc/html/index.html ./Documentation.html
