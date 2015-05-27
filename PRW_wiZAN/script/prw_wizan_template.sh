#!/bin/bash

#PBS -q production_qdr
#PBS -N prwizan
#PBS -l select=1:ncpus=1
#PBS -l place=free
#PBS -V

# change to the working directory
cd $PBS_O_WORKDIR
echo ">>>> Begin prw_wizan combination test"

#input file suffix
DIR='/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_maxTc/'
SUFFIX='N1Tc0.49to0.55'
EXT='.csv'
#file paths must be absolute pathways!!!!
FEAT='/scratch/hansaim.lim/wiZAN/PRW_wiZAN/input/zinc_chemIndex_feature.tsv'
TRAIN="${DIR}train_${SUFFIX}${EXT}"
TEST="${DIR}test_${SUFFIX}${EXT}"
PRW_OUT="${DIR}prw_out_temp_${SUFFIX}${EXT}"
PRWonly_OUT="${DIR}PRWonly_${SUFFIX}${EXT}"
TPR_OUT="${DIR}prwizan_TPR_${SUFFIX}${EXT}"

#parameters for PRW-wiZAN
PARA="[0.1, 300, 100, 0.75, 0.1]"

#parameters for wiZAN_dual.
#PARA="[0.001, 0.1, 0.01, 300, 100, 0.75, 0.5]"

# actual binary (with IO redirections) and required input
# parameters are called in the next line

#run PRW only test
python /scratch/hansaim.lim/wiZAN/PRW/PRW.py $TRAIN $TEST $FEAT $PRWonly_OUT 0.9

#run python PRW first
python /scratch/hansaim.lim/wiZAN/PRW_wiZAN/script/PRW_for_wiZAN.py $TRAIN $FEAT $PRW_OUT 0.9
#run matlab wiZAN with the output from PRW
matlab -r "cd /scratch/hansaim.lim/wiZAN/PRW_wiZAN/script/; PRW_wiZAN_onetest('$TRAIN', '$TEST', '$PRW_OUT','$TPR_OUT', '$PARA')"
#remove PRW output
rm $PRW_OUT
echo "Done!"
