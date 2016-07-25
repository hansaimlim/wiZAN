#!/bin/bash
#sleep 43200 #run after 12 hours

#
# should divide ikey files and run the pubchem search
#
#python ./chem_info_by_ikeys_pubchem.py ../x00 ../output/chem_info_x00.tsv > ../output/chem_x00_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x01 ../output/chem_info_x01.tsv > ../output/chem_x01_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x02 ../output/chem_info_x02.tsv > ../output/chem_x02_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x03 ../output/chem_info_x03.tsv > ../output/chem_x03_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x04 ../output/chem_info_x04.tsv > ../output/chem_x04_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x05 ../output/chem_info_x05.tsv > ../output/chem_x05_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x06 ../output/chem_info_x06.tsv > ../output/chem_x06_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x07 ../output/chem_info_x07.tsv > ../output/chem_x07_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x08 ../output/chem_info_x08.tsv > ../output/chem_x08_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x09 ../output/chem_info_x09.tsv > ../output/chem_x09_error.ext
#sleep 1200
#python ./chem_info_by_ikeys_pubchem.py ../x10 ../output/chem_info_x10.tsv > ../output/chem_x10_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x11 ../output/chem_info_x11.tsv > ../output/chem_x11_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x12 ../output/chem_info_x12.tsv > ../output/chem_x12_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x13 ../output/chem_info_x13.tsv > ../output/chem_x13_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x14 ../output/chem_info_x14.tsv > ../output/chem_x14_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x15 ../output/chem_info_x15.tsv > ../output/chem_x15_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x16 ../output/chem_info_x16.tsv > ../output/chem_x16_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x17 ../output/chem_info_x17.tsv > ../output/chem_x17_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x18 ../output/chem_info_x18.tsv > ../output/chem_x18_error.ext &
#python ./chem_info_by_ikeys_pubchem.py ../x19 ../output/chem_info_x19.tsv > ../output/chem_x19_error.ext

#sleep 1200 #wait 20 mins after the last job has done
#merge the files
#cat ../output/chem_info_x*.tsv >> ../output/chem_info_from_PubChem.tsv

python ./check_chems.py > ../output/check_chems_error.txt
python ./insert_chemicals.py > ../output/insert_chemicals_error.txt
python ./get_chem_disease_index.py > ../output/chem_disease_index_error.txt
python ./get_chem_prot_index.py > ../output/chem_prot_index_error.txt
