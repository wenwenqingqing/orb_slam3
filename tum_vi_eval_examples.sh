#!/bin/bash
pathDatasetTUM_VI='/Users/yyqing/Soft/DataSet/tum' #Example, it is necesary to change it by the dataset path


# Single Session Example

echo "------------------------------------"
echo "Evaluation of Magistrale 1 trajectory with Stereo-Inertial sensor"
python2 evaluation/evaluate_ate_scale.py "$pathDatasetTUM_VI"/dataset-corridor1_512_16/mav0/mocap0/data.csv f_dataset-corridor1_512_monoi.txt --plot corridor1_512_monoi.pdf
