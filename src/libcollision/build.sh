#!/bin/bash

g++ -c main.cpp 
g++ -c ugrid.cpp 
g++ ugrid.o main.o -o app

