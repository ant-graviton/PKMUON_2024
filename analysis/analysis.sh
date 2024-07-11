#!/bin/bash

NPROC=$(grep MHz /proc/cpuinfo | wc -l)
N=0

for ROOTFILE in $(ls -v ../build/root_file/CryMu_*.root); do
    [ $N = $NPROC ] && wait -n || let N+=1
    echo root -l -q -b 'src/analysis.cc("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAna}"'")'
    root -l -q -b 'src/analysis.cc("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAna}"'")' &>${ROOFILE/CryMu/CryMuAna}.log &
done
wait
hadd -f ../build/root_file/CryMuAna.root ../build/root_file/CryMuAna_*.root