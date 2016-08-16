#!/bin/bash

#PBS -q production_qdr
#PBS -N NRLMF_N1Tc
#PBS -l select=2:ncpus=1:mem=8000mb
#PBS -l place=free
#PBS -V

# change to the working directory
cd $PBS_O_WORKDIR
datadir_pre="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/"
echo ">>>> Begin NRLMF"
#python PyDTI_10cv.py --method="nrlmf" --dataset="N1Tc0.5to0.6" --data-dir="${datadir_pre}NTMaxTc/N1/Tc5_6/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N1Tc0.5to0.6_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N1Tc0.6to0.7" --data-dir="${datadir_pre}NTMaxTc/N1/Tc6_7/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N1Tc0.6to0.7_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N1Tc0.7to0.8" --data-dir="${datadir_pre}NTMaxTc/N1/Tc7_8/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N1Tc0.7to0.8_r100i300.txt

#python PyDTI_10cv.py --method="nrlmf" --dataset="N1Tc0.8to0.9" --data-dir="${datadir_pre}NTMaxTc/N1/Tc8_9/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N1Tc0.8to0.9_r100i300.txt

python PyDTI_10cv.py --method="nrlmf" --dataset="N1Tc0.9to1.0" --data-dir="${datadir_pre}NTMaxTc/N1/Tc9_10/" --cvs=1 --predict-num=5 --method-opt="r=100 max_iter=300" >> /scratch/hansaim.lim/wiZAN/PyDTI/NRLMF_on_ZINC/10cv/N1Tc0.9to1.0_r100i300.txt
