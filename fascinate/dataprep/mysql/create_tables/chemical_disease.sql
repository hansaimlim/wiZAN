CREATE TABLE chemical_disease(
chemical_index INT(10) UNSIGNED  NOT NULL,
disease_index INT(10) UNSIGNED NOT NULL,
PRIMARY KEY(chemical_index, disease_index),
UNIQUE KEY(chemical_index, disease_index));
