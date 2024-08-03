#!/bin/bash

NPROC=$(nproc || sysctl -n hw.logicalcpu || getconf _NPROCESSORS_ONLN)
N=0
PIDS=()

for ROOTFILE in $(ls -v ../build/root_file/CryMu_*.root); do
    if [ $N = $NPROC ]; then
        wait $PIDS
        PIDS=(${PIDS[@]:1})
    else
        let N+=1
    fi
    (
        echo root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAna}"'")'
        root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAna}"'")' &>${ROOTFILE/CryMu/CryMuAna}.log
        echo root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAll}"'", true)'
        root -l -q -b 'src/analysis.C("'"${ROOTFILE}"'", "'"${ROOTFILE/CryMu/CryMuAll}"'", true)' &>${ROOTFILE/CryMu/CryMuAll}.log
    ) &
    PIDS=(${PIDS[@]} $!)
done
wait
(
    hadd -f ../build/root_file/CryMuAna.root $(ls -v ../build/root_file/CryMuAna_*.root)
    root -l -q -b 'src/PoCA_sim.C("../build/root_file/CryMuAna.root", "../build/root_file/CryMuAnaPoCA.root")'
    root -l -q -b 'src/draw_sim.C("../build/root_file/CryMuAnaPoCA.root", "draw_sim_ana.pdf")'
) &
(
    hadd -f ../build/root_file/CryMuAll.root $(ls -v ../build/root_file/CryMuAll_*.root)
    root -l -q -b 'src/PoCA_sim.C("../build/root_file/CryMuAll.root", "../build/root_file/CryMuAllPoCA.root")'
    root -l -q -b 'src/draw_sim.C("../build/root_file/CryMuAllPoCA.root", "draw_sim_all.pdf")'
) &
wait
