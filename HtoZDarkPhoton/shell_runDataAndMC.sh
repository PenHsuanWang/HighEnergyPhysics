#!/bin/bash 

nohup root -l -b -q runData.C &> test_runMC_Da.txt &
nohup root -l -b -q runMC_background.C &> test_runMC_BG.txt &
