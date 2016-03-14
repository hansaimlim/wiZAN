#!/usr/bin/python
# Comparative modeling with multiple templates
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '/home/hlim/wiZAN/RhoGAP/templates/']

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Add some restraints from a file:
#       rsr.append(file='my_rsrs1.rsr')

#       Residues 20 through 30 should be an alpha helix:
        rsr.add(secondary_structure.alpha(self.residue_range('25:', '31:')))
        rsr.add(secondary_structure.alpha(self.residue_range('41:', '53:')))
        rsr.add(secondary_structure.alpha(self.residue_range('67:', '77:')))
        rsr.add(secondary_structure.alpha(self.residue_range('89:', '102:')))
        rsr.add(secondary_structure.alpha(self.residue_range('110:', '121:')))
        rsr.add(secondary_structure.alpha(self.residue_range('125:', '138:')))
        rsr.add(secondary_structure.alpha(self.residue_range('142:', '160:')))
        rsr.add(secondary_structure.alpha(self.residue_range('172:', '175:')))
        rsr.add(secondary_structure.alpha(self.residue_range('188:', '210:')))
#       Two beta-strands:
        rsr.add(secondary_structure.strand(self.residue_range('17:', '22:')))
        rsr.add(secondary_structure.strand(self.residue_range('59:', '63:')))
#       An anti-parallel sheet composed of the two strands:
#        rsr.add(secondary_structure.sheet(at['N:1'], at['O:14'],
#                                          sheet_h_bonds=-5))
#       Use the following instead for a *parallel* sheet:
#       rsr.add(secondary_structure.sheet(at['N:1'], at['O:9'],
#                                         sheet_h_bonds=5))

#       Restrain the specified CA-CA distance to 10 angstroms (st. dev.=0.1)
#       Use a harmonic potential and X-Y distance group.
#        rsr.add(forms.gaussian(group=physical.xy_distance,
#                               feature=features.distance(at['CA:35'],
#                                                         at['CA:40']),
#                               mean=10.0, stdev=0.1))

a = MyModel(env,
              alnfile  = '/home/hlim/wiZAN/RhoGAP/comp155_c0_seq1/align/comp155_c0_seq1_noSig.pir', # alignment filename
              knowns   = ('2mbg','2ovj','3wpq'),     # codes of the templates
							#template 2mbg from I-TASSER; 2ovj and 3wpq from Lb reference
              sequence = 'comp155_c0_seq1')               # code of the target
a.starting_model= 701                 # index of the first model
a.ending_model  = 705               # index of the last model
                                    # (determines how many models to calculate)
#aln = alignment(env)
#aln.append(file='clustalo-Lbtemp_pdbseq.pir', align_codes='all')
#aln.check()

#Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 500

#thorough MD optimization:
a.md_level = refine.slow

#repeat the whole cycle 5 times and do not stop unless obj.func. > 1e6
a.repeat_optimization = 100
a.max_molpdf = 1e6

a.make()                            # do the actual comparative modeling
