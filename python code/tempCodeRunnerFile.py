# %%
from Bio.PDB import *
import nglview 
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('2FAT', '2FAT.cif')
view = nglview.show_biopython(structure)
view.clear_representations()
#view as ball and stick (atom and bond)
view.add_ball_and_stick()
view