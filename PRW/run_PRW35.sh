#!/bin/bash
x=1
while [ $x -le 20 ]
do
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N1L1to5.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N1L1to5.csv ./tpr35/TPR35_N1L1to5.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N1L6to10.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N1L6to10.csv ./tpr35/TPR35_N1L6to10.txt 0.9 & 
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N1L11to15.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N1L11to15.csv ./tpr35/TPR35_N1L11to15.txt 0.9 & 
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N1L16to20.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N1L16to20.csv ./tpr35/TPR35_N1L16to20.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N1L21more.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N1L21more.csv ./tpr35/TPR35_N1L21more.txt 0.9 &

	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N1Tc0.49to0.6.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N1Tc0.49to0.6.csv ./tpr35/TPR35_N1Tc0.49to0.6.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N1Tc0.6to0.7.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N1Tc0.6to0.7.csv ./tpr35/TPR35_N1Tc0.6to0.7.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N1Tc0.7to0.8.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N1Tc0.7to0.8.csv ./tpr35/TPR35_N1Tc0.7to0.8.txt 0.9 &

	wait
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N1Tc0.8to0.9.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N1Tc0.8to0.9.csv ./tpr35/TPR35_N1Tc0.8to0.9.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N1Tc0.9to1.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N1Tc0.9to1.csv ./tpr35/TPR35_N1Tc0.9to1.txt 0.9 &


	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N2L1to5.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N2L1to5.csv ./tpr35/TPR35_N2L1to5.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N2L6to10.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N2L6to10.csv ./tpr35/TPR35_N2L6to10.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N2L11to15.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N2L11to15.csv ./tpr35/TPR35_N2L11to15.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N2L16to20.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N2L16to20.csv ./tpr35/TPR35_N2L16to20.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N2L21more.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N2L21more.csv ./tpr35/TPR35_N2L21more.txt 0.9 &

	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N2Tc0.49to0.6.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N2Tc0.49to0.6.csv ./tpr35/TPR35_N2Tc0.49to0.6.txt 0.9 &
	wait

	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N2Tc0.6to0.7.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N2Tc0.6to0.7.csv ./tpr35/TPR35_N2Tc0.6to0.7.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N2Tc0.7to0.8.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N2Tc0.7to0.8.csv ./tpr35/TPR35_N2Tc0.7to0.8.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N2Tc0.8to0.9.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N2Tc0.8to0.9.csv ./tpr35/TPR35_N2Tc0.8to0.9.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N2Tc0.9to1.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N2Tc0.9to1.csv ./tpr35/TPR35_N2Tc0.9to1.txt 0.9 &

	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N3upL1to5.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N3upL1to5.csv ./tpr35/TPR35_N3upL1to5.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N3upL6to10.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N3upL6to10.csv ./tpr35/TPR35_N3upL6to10.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N3upL11to15.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N3upL11to15.csv ./tpr35/TPR35_N3upL11to15.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N3upL16to20.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N3upL16to20.csv ./tpr35/TPR35_N3upL16to20.txt 0.9 &

	wait
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_numligand/train_N3upL21more.csv ../ZINC_data/chem_prot/numtarget_numligand/test_N3upL21more.csv ./tpr35/TPR35_N3upL21more.txt 0.9 &

	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N3upTc0.49to0.6.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N3upTc0.49to0.6.csv ./tpr35/TPR35_N3upTc0.49to0.6.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N3upTc0.6to0.7.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N3upTc0.6to0.7.csv ./tpr35/TPR35_N3upTc0.6to0.7.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N3upTc0.7to0.8.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N3upTc0.7to0.8.csv ./tpr35/TPR35_N3upTc0.7to0.8.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N3upTc0.8to0.9.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N3upTc0.8to0.9.csv ./tpr35/TPR35_N3upTc0.8to0.9.txt 0.9 &
	python PRW_TPR35.py ../ZINC_data/chem_prot/numtarget_maxTc/widerange/train_N3upTc0.9to1.csv ../ZINC_data/chem_prot/numtarget_maxTc/widerange/test_N3upTc0.9to1.csv ./tpr35/TPR35_N3upTc0.9to1.txt 0.9 &
	wait
	x=$(( $x + 1 ))
done
