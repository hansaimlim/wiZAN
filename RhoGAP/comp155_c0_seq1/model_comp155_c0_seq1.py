#!/usr/bin/python
# Comparative modeling with multiple templates
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '/home/hansaim/wiZAN/RhoGAP/templates/']

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
#       Add some restraints from a file:
#       rsr.append(file='my_rsrs1.rsr')

#       secondary structures predicted by majority consensus of SYMPRED, SSPRO, PsiPred, YASPIN, JPred, numbering without signaling peptide region
        rsr.add(secondary_structure.alpha(self.residue_range('6:', '17:')))
        rsr.add(secondary_structure.alpha(self.residue_range('33:', '42:')))
        rsr.add(secondary_structure.alpha(self.residue_range('54:', '66:')))
        rsr.add(secondary_structure.alpha(self.residue_range('78:', '85:')))
        rsr.add(secondary_structure.alpha(self.residue_range('90:', '103:')))
        rsr.add(secondary_structure.alpha(self.residue_range('109:', '124:')))
        rsr.add(secondary_structure.alpha(self.residue_range('153:', '175:')))

#       beta-strands:
        rsr.add(secondary_structure.strand(self.residue_range('137:', '140:')))

a = MyModel(env,
              alnfile  = '/home/hansaim/wiZAN/RhoGAP/comp155_c0_seq1/comp155_c0_seq1.pir', # alignment filename
              knowns   = ('3wpq','3byi','3iug'),     
              sequence = 'comp155_c0_seq1')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 3               # index of the last model
                                    # (determines how many models to calculate)
#Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 500

#thorough MD optimization:
a.md_level = refine.slow

#repeat the whole cycle 10 times and do not stop unless obj.func. > 1e6
a.repeat_optimization = 10
a.max_molpdf = 1e6

a.make()                            # do the actual comparative modeling
