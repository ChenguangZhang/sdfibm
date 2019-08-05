#!/bin/bash
echo " ======= make shape ======= "
cd libshape
make clean
make 
cd ..

echo " ======= make motion ======= "
cd libmotion
make clean
make 
cd ..

echo " ======= make collision ======= "
cd libcollision
make clean
make 
cd ..

echo " ======= make sdfibm ======= "
make clean
make 

