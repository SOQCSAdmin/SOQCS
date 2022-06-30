#!/bin/bash

##### CONFIGURATION ########
VER=3.4.0     #VERSION
DIR=./cache   #CACHE DIRECTORY
INSTALL=./src/Eigen  #EIGEN INSTALL DIRECTION (DON'T TOUCH)

##### SCRIPT ########
mkdir $DIR
wget -O $DIR/eigen.tar.gz  https://gitlab.com/libeigen/eigen/-/archive/$VER/eigen-$VER.tar.gz
tar xvf $DIR/eigen.tar.gz -C $DIR
cp -R $DIR/eigen-$VER/Eigen/* $INSTALL/
ln -s ./src ./py_src/src
