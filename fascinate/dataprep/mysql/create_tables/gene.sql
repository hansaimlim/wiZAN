CREATE TABLE gene(
gene_index INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
HGNC INT(10) UNSIGNED NOT NULL,
gene_symbol VARCHAR(50),
uniprot_accession VARCHAR(30),
ensembl VARCHAR(50),
protein_name VARCHAR(255),
chromosome VARCHAR(50),
sequence VARCHAR(40000),
sequence_source VARCHAR(30),
PRIMARY KEY(gene_index),
UNIQUE KEY(HGNC));
