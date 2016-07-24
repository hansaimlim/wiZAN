CREATE TABLE gene_disease(
gene_index INT(10) UNSIGNED NOT NULL,
disease_index INT(10) UNSIGNED NOT NULL,
status VARCHAR(255),
PRIMARY KEY(gene_index, disease_index),
UNIQUE KEY(gene_index, disease_index));
