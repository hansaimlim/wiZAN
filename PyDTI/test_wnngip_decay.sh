#!/bin/sh

for decay in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 
do
	python PyDTI.py --method="wnngip" --dataset="N1L1to5" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  > /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N1L6to10" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N1L11to15" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N1L16to20" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N1L21more" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt

	python PyDTI.py --method="wnngip" --dataset="N2L1to5" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N2L6to10" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N2L11to15" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N2L16to20" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N2L21more" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt

	python PyDTI.py --method="wnngip" --dataset="N3upL1to5" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N3upL6to10" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N3upL11to15" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N3upL16to20" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
	python PyDTI.py --method="wnngip" --dataset="N3upL21more" --data-dir="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/" --cvs=1 --predict-num=5 --method-opt="T=$decay"  >> /scratch/hansaim.lim/PyDTI/NRLMF_on_ZINC/wnngip_decay.txt
done
