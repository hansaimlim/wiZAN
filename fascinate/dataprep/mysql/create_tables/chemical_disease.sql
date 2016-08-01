CREATE TABLE chemical_disease(
chemical_disease_index INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
chemical_index INT(10) UNSIGNED  NOT NULL,
disease_index INT(10) UNSIGNED NOT NULL,
PRIMARY KEY(chemical_disease_index),
UNIQUE KEY(chemical_index, disease_index));
