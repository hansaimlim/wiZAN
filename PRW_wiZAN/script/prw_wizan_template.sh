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
DirIN="/scratch/hansaim.lim/wiZAN/ZINC_data/chem_prot/numtarget_numligand/"
DirOUT="/scratch/hansaim.lim/prw_wizan_combination/numtarget_numligand/"
SUFFIX="N1L1to5"
EXT=".csv"
EXT_out=".tsv"
#file paths must be absolute pathways!!!!
FEAT="/scratch/hansaim.lim/wiZAN/PRW_wiZAN/input/zinc_chemIndex_feature.tsv"
TRAIN=${DirIN}"train_"${SUFFIX}${EXT}
TEST=${DirIN}"test_"${SUFFIX}${EXT}
PRW_OUT=${DirOUT}"prw_out_temp_"${SUFFIX}${EXT}	#csv output
PRWonly_OUT=${DirOUT}"PRWonly_"${SUFFIX}${EXT_out}	#non csv
wiZAN_OUT=${DirOUT}"wiZAN_dual_"${SUFFIX}${EXT_out}	#non csv
PRWIZAN_OUT=${DirOUT}"prwizan_TPR_"${SUFFIX}${EXT_out}	#non csv

# actual binary (with IO redirections) and required input
# parameters are called in the next line

#run PRW only test
python /scratch/hansaim.lim/wiZAN/PRW/PRW.py $TRAIN $TEST $FEAT $PRWonly_OUT 0.9

#run python PRW first
python /scratch/hansaim.lim/wiZAN/PRW_wiZAN/script/PRW_for_wiZAN.py $TRAIN $FEAT $PRW_OUT 0.9
#run matlab wiZAN with the output from PRW
matlab -r "cd /scratch/hansaim.lim/wiZAN/PRW_wiZAN/script/; PRW_wiZAN_onetest('$TRAIN', '$TEST', '$PRW_OUT', '$PRWIZAN_OUT')"
#remove PRW output
rm $PRW_OUT
#run wiZAN_dual
matlab -r "cd /scratch/hansaim.lim/wiZAN/wiZAN_dual/; wiZAN_dual_csv('$TRAIN', '$TEST', '$wiZAN_OUT')"
echo "Done!"
