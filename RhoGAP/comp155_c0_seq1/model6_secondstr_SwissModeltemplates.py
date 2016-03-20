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

#       secondary structures predicted by SYMPRED, without signaling peptide region
        rsr.add(secondary_structure.alpha(self.residue_range('4:', '10:')))
        rsr.add(secondary_structure.alpha(self.residue_range('20:', '32:')))
        rsr.add(secondary_structure.alpha(self.residue_range('46:', '56:')))
        rsr.add(secondary_structure.alpha(self.residue_range('68:', '81:')))
        rsr.add(secondary_structure.alpha(self.residue_range('89:', '100:')))
        rsr.add(secondary_structure.alpha(self.residue_range('104:', '117:')))
        rsr.add(secondary_structure.alpha(self.residue_range('121:', '139:')))
        rsr.add(secondary_structure.alpha(self.residue_range('167:', '189:')))
        rsr.add(secondary_structure.alpha(self.residue_range('242:', '244:')))
#       Two beta-strands:
        rsr.add(secondary_structure.strand(self.residue_range('38:', '42:')))
        rsr.add(secondary_structure.strand(self.residue_range('151:', '154:')))

a = MyModel(env,
              alnfile  = '/home/hlim/wiZAN/RhoGAP/comp155_c0_seq1/align/comp155_c0_seq1_SwissModel_align.pir', # alignment filename
              knowns   = ('3wpq','3wps','3w6r','3cxl','2mbg'),     # first 4 templates high scored by SwissModel, 5th for c-term
							#template 2mbg from I-TASSER;
              sequence = 'comp155_c0_seq1')               # code of the target
a.starting_model= 801                 # index of the first model
a.ending_model  = 805               # index of the last model
                                    # (determines how many models to calculate)
#Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 500

#thorough MD optimization:
a.md_level = refine.slow

#repeat the whole cycle 100 times and do not stop unless obj.func. > 1e6
a.repeat_optimization = 100
a.max_molpdf = 1e6

a.make()                            # do the actual comparative modeling
