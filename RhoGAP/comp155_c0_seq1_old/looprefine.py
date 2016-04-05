# Loop refinement of an existing model
from modeller import *
from modeller.automodel import *
#from modeller import soap_loop

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '/home/hansaim/RhoGAP/templates/']
# Create a new class based on 'loopmodel' so that we can redefine
# select_loop_atoms (necessary)
class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # Three loops simultaneously
        return selection(self.residue_range('13:', '34:'),
                         self.residue_range('166:', '178:'),
                         self.residue_range('199:','217:'))

m = MyLoop(env,
           inimodel='/home/hansaim/RhoGAP/comp155_c0_seq1/align/comp155_c0_seq1.B99990001.pdb',   # initial model of the target
           sequence='comp155_c0_seq1',                 # code of the target
           loop_assess_methods=assess.DOPE) # assess loops with DOPE
#          loop_assess_methods=soap_loop.Scorer()) # assess with SOAP-Loop

m.loop.starting_model= 22           # index of the first loop model
m.loop.ending_model  = 23           # index of the last loop model

#Very thorough VTFM optimization:
m.library_schedule = autosched.slow
m.max_var_iterations = 500
m.loop.md_level = refine.slow  # loop refinement method
#repeat the whole cycle 2 times and do not stop unless obj.func. > 1e6
m.repeat_optimization = 3
m.max_molpdf = 1e6

m.make()
