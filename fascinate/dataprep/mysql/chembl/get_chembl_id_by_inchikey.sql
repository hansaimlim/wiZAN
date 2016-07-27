SELECT mdic.chembl_id, mdic.pref_name, cstr.standard_inchi_key, cstr.canonical_smiles
FROM chembl_20.compound_structures cstr
INNER JOIN chembl_20.molecule_dictionary mdic ON cstr=molregno=mdic=molregno
WHERE cstr.standard_inchi_key=%s
