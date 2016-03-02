#!/usr/bin/python
# Comparative modeling with multiple templates
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '/home/hansaim/RhoGAP/templates/']

a = automodel(env,
              alnfile  = './align/clustalo-Lbtemp_pdbseq.pir.bak', # alignment filename
              knowns   = ('1xa6', '1f7c','1he1'),     # codes of the templates
              sequence = 'comp155_c0_seq1')               # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)
#aln = alignment(env)
#aln.append(file='clustalo-Lbtemp_pdbseq.pir', align_codes='all')
#aln.check()

#Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

#thorough MD optimization:
a.md_level = refine.slow

#repeat the whole cycle 2 times and do not stop unless obj.func. > 1e6
a.repeat_optimization = 3
a.max_molpdf = 1e6

a.make()                            # do the actual comparative modeling
