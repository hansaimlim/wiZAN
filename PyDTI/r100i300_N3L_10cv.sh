#!/bin/bash

#PBS -q production_qdr
#PBS -N NRLMF_N3L
#PBS -l select=2:ncpus=1:mem=8000mb
#PBS -l place=free
#PBS -V

# change to the working directory
cd $PBS_O_WORKDIR
datadir_pre="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/"
echo ">>>> Begin NRLMF"
#python PyDTI_10cv.py --method="nrlmf" --dataset="N3L1to5" --data-dir="${datadir_pre}NTNL/N3/L1to5/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N3L1to5_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N3L6to10" --data-dir="${datadir_pre}NTNL/N3/L6to10/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N3L6to10_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N3L11to15" --data-dir="${datadir_pre}NTNL/N3/L11to15/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N3L11to15_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N3L16to20" --data-dir="${datadir_pre}NTNL/N3/L16to20/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N3L16to20_r100i300.txt

python PyDTI_10cv.py --method="nrlmf" --dataset="N3L21more" --data-dir="${datadir_pre}NTNL/N3/L21more/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N3L21more_r100i300.txt
