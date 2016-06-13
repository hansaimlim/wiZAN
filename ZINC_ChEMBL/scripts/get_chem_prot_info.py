#!/usr/bin/python

import MySQLdb as db
import sys

con = db.connect('localhost', 'root', 'rhkdlf2043', 'chembl_20');
cur=con.cursor()

infile='./ZC_ambiguous_pairs.tsv'
print "Standard_InChI_Key\tUniProt\tIC50\tUnit\tType"
for line in open(infile,"r").xreadlines():
        z=line.strip().split("\t")
        ikey=str(z[0])
        uniprot=str(z[1])

        qry="SELECT comp_str.standard_inchi_key, com_seq.accession, ac.standard_value, ac.standard_units, ac.standard_type\
             FROM compound_structures comp_str\
	     INNER JOIN activities ac ON comp_str.molregno=ac.molregno\
             INNER JOIN assays ay ON ac.assay_id=ay.assay_id\
             INNER JOIN target_components tc ON ay.tid=tc.tid\
             INNER JOIN component_sequences com_seq ON tc.component_id=com_seq.component_id\
             WHERE ac.standard_type='IC50' AND comp_str.standard_inchi_key='%s'" % (ikey)
        cur.execute(qry)
        result=cur.fetchall()
	for r in result:
		ikey=str(r[0])
		unip=str(r[1])
		conc=str(r[2])
		unit=str(r[3])
		assay=str(r[4])
		print "%s\t%s\t%s\t%s\t%s" % (ikey,unip,conc,unit,assay)

con.close()
