CREATE TABLE CAS_InChIKey(
CAS_index INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
CAS VARCHAR(40) NOT NULL,
InChIKey VARCHAR(40) NOT NULL,
PRIMARY KEY(CAS_index),
UNIQUE KEY(CAS, InChIKey));