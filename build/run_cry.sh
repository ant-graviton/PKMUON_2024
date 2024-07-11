#!/bin/bash

mkdir -p root_file

# 循环运行Geant4模拟，传入不同的i值
for i in $(seq 1 $(lscpu | grep '^CPU(s):' | grep -o '[0-9]\+')); do

  cp CryMu.mac CryMu_$i.mac
  sed -i "s/CryMu\\.root/CryMu_${i}.root/g" CryMu_$i.mac

  # 运行Geant4模拟
  ./muPos CryMu_$i.mac &

done

#rm -rf CryMu_*.mac
