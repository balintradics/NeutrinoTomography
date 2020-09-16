#!/bin/sh
# lxplus batch script 

export OUTPATH="/afs/cern.ch/work/b/bradics/public/nutomo"
export SRCPATH="/afs/cern.ch/work/b/bradics/Spark/nutomo"
export LOCALPATH="$PWD"

echo $LOCALPATH

echo "Phase 1: build"
# copy executable to pool PWD and build
cp -r $SRCPATH .
mkdir build
cd build 
cmake3 ../nutomo
make

echo "Phase 2: run exe"
# Run the exe
cd bin
./Run_NuProb_OscLongitude

echo "Phase 3: mv data to output path"
cp ./NuProb_300km.txt $OUTPATH


