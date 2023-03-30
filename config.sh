#!/bin/bash

prompt=true
while [ $prompt = "true" ]; 
do
    read -p "Choose your hardware:
             [Dd] Dev/Manual (Linux with OpenNP) 
             [Ll] Linux
             [Mm] MacOsX (native)
             [Cc] MacOsx (conda)
              " hdw
             
    case $hdw in
        [Dd]* ) cp conf.inc_man conf.inc 
                cp download_man.sh download.sh
                prompt=false;;
        [Ll]* ) cp conf.inc_lin conf.inc 
                cp download_lin.sh download.sh
                prompt=false;;
        [Mm]* ) cp conf.inc_osx conf.inc 
                cp download_osx.sh download.sh
                prompt=false;;
        [Cc]* ) cp conf.inc_conda conf.inc 
                cp download_osx.sh download.sh
                prompt=false;;
        * )     echo "Invalid option" 
                echo " " 
    esac
done

chmod 744 *.sh
./download.sh
./links.sh
cd py_src
 chmod 744 *.sh
./py_conf.sh
cd ..
