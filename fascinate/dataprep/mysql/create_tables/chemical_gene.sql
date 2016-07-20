CREATE TABLE chemical_gene(
chemical_index INT(10) UNSIGNED NOT NULL,
gene_index INT(10) UNSIGNED NOT NULL,
source VARCHAR(100) NOT NULL,
PRIMARY KEY(chemical_index, gene_index),
UNIQUE KEY(chemical_index, gene_index));
