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
        rsr.add(secondary_structure.alpha(self.residue_range('9:', '20:')))
        rsr.add(secondary_structure.alpha(self.residue_range('35:', '45:')))
        rsr.add(secondary_structure.alpha(self.residue_range('57:', '69:')))
        rsr.add(secondary_structure.alpha(self.residue_range('81:', '88:')))
        rsr.add(secondary_structure.alpha(self.residue_range('93:', '106:')))
        rsr.add(secondary_structure.alpha(self.residue_range('112:', '127:')))
        rsr.add(secondary_structure.alpha(self.residue_range('155:', '179:')))

#       beta-strands:
        rsr.add(secondary_structure.strand(self.residue_range('140:', '143:')))

a = MyModel(env,
              alnfile  = '/home/hansaim/wiZAN/RhoGAP/comp155_c0_seq3/comp155_c0_seq3.pir', # alignment filename
              knowns   = ('3wps','3byi','3iug','1f7c'),     
              sequence = 'comp155_c0_seq3')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 5               # index of the last model
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
