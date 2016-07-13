The java code used to process raw data and reconstruct inferred dependency is given in this code. 
Its sample output for indexed layer network is in 'chem_chem_tg_sim0.5.tsvOut'.

The sample indexed dependency file is given in 'chem_gene.tsvdep', which is based on the name-index map given in 'chem_chem_tg_sim0.5.tsvMap' and 'gene_gene_ppi.tsvMap'.

To reconstruct the indexed dependency back to entity dependency, the format of the dependency file should be the same with the 'chem_gene.tsvdep'. (i.e., indices are seperated by '\t', each pair of dependency is seperated by '\n').