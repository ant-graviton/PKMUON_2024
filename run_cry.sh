#!/bin/bash

NPROC=$(grep MHz /proc/cpuinfo | wc -l)
mkdir -p root_file

# 循环运行Geant4模拟，传入不同的i值
for i0 in $(seq 0 11); do
  for i1 in $(seq 0 $[NPROC-1]); do
    i=$(($i0 * $NPROC + $i1))
    cp CryMu.mac CryMu_$i.mac
    sed -i "s/CryMu\\.root/CryMu_${i}.root/g" CryMu_$i.mac

    # 运行Geant4模拟
    ./muPos CryMu_$i.mac &
  done
  wait
done

#rm -rf CryMu_*.mac
