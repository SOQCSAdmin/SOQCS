#!/bin/bash

chmod 744 *.sh
./download.sh
./links.sh
cd py_src
 chmod 744 *.sh
./py_conf.sh
cd ..
