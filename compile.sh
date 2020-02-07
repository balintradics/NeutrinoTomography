#!/bin/bash

g++ -c neutrino_osc.cpp
g++ -c complex.cpp
g++ -c thematrix.cpp
g++ -c nrutil.cpp
g++ -c odeint.cpp
g++ -o testme testme.cpp neutrino_osc.o complex.o thematrix.o nrutil.o odeint.o
