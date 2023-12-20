#!/bin/bash
cd $CMSSW_BASE/src/CSCLayerAlignment/CSCLayerAlignmentProducer/scripts/
eval `scramv1 runtime -sh`
g++ -g -o CSCLayerAlignmentFitter.o CSCLayerAlignmentFitter.cpp $(root-config --cflags --glibs --ldflags) -lMinuit
./CSCLayerAlignmentFitter.o
