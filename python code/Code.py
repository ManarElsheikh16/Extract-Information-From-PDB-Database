from Bio.PDB import *

#local
file=open("pdb2fat.ent","r")
header_dict= parse_pdb_header(file)
print(list(header_dict))
print(header_dict)
file.close()

#from server
pdb=PDBList()
pdb.retrieve_pdb_file("2FAT",pdir=".",file_format="pdb")
parser=PDBParser(PERMISSIVE=True,QUIET=True) #the PERMISSIVE flag indicates that-- a number of common problems associated with PDB files will be ignored
data=parser.get_structure("2FAT","pdb2FAT.ent")

print(data)                               #print id of file
print(type(data))                         #print type of data
print(data.header.keys())                 #print headers of file
print(data.header["deposition_date"])     #print the content of deposition_date
print(data.header["release_date"])        #print the content of release_date
print(data.header["name"])                #print the content of name
print(data.header["head"])                #print the content of head
print(data.header["idcode"])              #print the content of idcode
print(data.header["structure_method"])    #print the content of structure_method
print(data.header["resolution"])          #print the content of resolution
print(data.header["structure_reference"]) #print the content of structure_reference
print(data.header["journal_reference"])   #print the content of journal_reference
print(data.header["author"])              #print the content of author
print(data.header["compound"])            #print the content of compound
print(data.header["source"])              #print the content of source
print(data.header["keywords"])            #print the content of keywords
print(data.header["journal"])             #print the content of journal
print(data.header["missing_residues"])    #print the content of missing_residues
print(data.header["has_missing_residues"])#print the content of has_missing_residues

#get the structure of the atom
model=data.get_models()
print(model)                               #print model of data
models=list(model)
print(models)                              #print id of model
print(type(models[0]))                     #print type of model

chain=models[0].get_chains()
chains=list(chain)                         
print(chains)                              #print chains exist in model 
print(type(chains[0]))                     #print type of chains

residue_chain_L=chains[0].get_residues()   
residues_chain_L=list(residue_chain_L)
print(len(residues_chain_L))               #print the length of residues exist in chain_L
print(residues_chain_L)                    #print residues exist in chain_L

residue_chain_H=chains[1].get_residues()
residues_chain_H=list(residue_chain_H)    
print(len(residues_chain_H))              #print the length of residues exist in chain_H
print(residues_chain_H)                   #print residues exist in chain_H

atom=residues_chain_L[0].get_atoms()
atoms=list(atom)
print(atoms)                              #print the list of atoms exist in residues_chain_L
print(atoms[7].get_vector())              #print the x,y,z coordinate of atom in index (7)

atom=residues_chain_H[1].get_atoms()
atoms=list(atom)
print(atoms)                              #print the list of atoms exist in residues_chain_H
print(atoms[6].get_vector())              #print the x,y,z coordinate of atom in index (6)

###########################################################################################
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()
structure = parser.get_structure("2FAT", "pdb2FAT.ent")
model = structure[0]
chain=model.get_chains()
###########################################################################################
from Bio.PDB.PDBParser import PDBParser   #1
p = PDBParser()
structure = p.get_structure('X', 'pdb2fat.ent')
for model in structure:
     for chain in model:
         for residue in chain:
             for atom in residue:
                 print(atom)                   #3 methods to print all atoms exist in structure               
###########################################################################################                
from Bio.PDB.PDBParser import PDBParser  #2
p = PDBParser()
structure = p.get_structure('X', 'pdb2fat.ent')
atoms = structure.get_atoms()
for atom in atoms:
     print(atom)
###########################################################################################    
atoms = chain.get_atoms()  #3 
for atom in atoms:
     print(atom)
     ####################################### wait
#get all residues from a structure: A=atom, R=residue, C=chain, M=model, S=structure
chain_list = Selection.unfold_entities(model, 'C')
print(chain_list)                                  #print id of chains that exist in model
atom_list = Selection.unfold_entities(chain, 'A')
print(len(atom_list))                              #print length of atoms found in chain
print(atom_list)                                   #print all atoms in chains
###########################################################################################
residues = model.get_residues()  #1
for residue in residues:
     print(residue)                                     #print all residues exist in chain L and chain H

res_list = Selection.unfold_entities(structure, 'R')  #2
print(res_list)
###########################################################################################
for model in structure.get_list():
     for chain in model.get_list():
         for residue in chain.get_list():
             if residue.has_id("CA"):
                 ca = residue["CA"]
                 if ca.get_bfactor() > 45.0:
                     print(ca.get_coord())    #print coordinate of x,y,z to the residues that has id (CA) and bfactor > 45.0
###########################################################################################
for model in structure.get_list():
     for chain in model.get_list():
         for residue in chain.get_list():
             if residue.is_disordered():
                 res_id = residue.get_id()[1]
                 resname = residue.get_resname()
                 model_id = model.get_id()
                 chain_id = chain.get_id()
                 print(model_id, chain_id, resname, res_id) #print disorder residues
else:
    print("all residues does not contain disordered atoms")
###########################################################################################
for model in structure.get_list():
     for chain in model.get_list():
         for residue in chain.get_list():
            for atom in residue.get_list():
                if atom.is_disordered():
                    print(atom)                #print disorder atom
else:
    print("all atoms is not disordered")
###########################################################################################
for residue in chain.get_list():
    residue_id = residue.get_id()
    hetfield = residue_id[0]
    if hetfield[0]=="W":
        print(residue_id)                #Print all water residues in chains
################################## mmcif file #####################################
from Bio.PDB import *
pdb_2=PDBList()
pdb_2.retrieve_pdb_file('2FAT' , pdir = '.', file_format = 'mmCif')
parser_2=MMCIFParser(QUIET=True)
data_2=parser_2.get_structure('2FAT','2FAT.cif')

print(data_2)                              #print the structure id of the file
print(type(data_2))                        #print type of data
print(data_2.header.keys())                #print headers of file
print(data_2.header["deposition_date"])    #print the content of deposition_data   
print(data_2.header["name"])               #print the content of name
print(data_2.header["head"])               #print the content of head
print(data_2.header["idcode"])             #print the content of idcode
print(data_2.header["structure_method"])   #print the content of structure_method
print(data_2.header["resolution"])         #print the content of resolution
###########################################################################################
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
mmcif_dict=MMCIF2Dict('2FAT.cif')
#print sheet_id,range_id_1,range_id_2,label_atom_id and label_comp_id of hydrogen bond from the backbone of atoms and amino acids
sc_1=mmcif_dict["_pdbx_struct_sheet_hbond.sheet_id"]
sc_2=mmcif_dict["_pdbx_struct_sheet_hbond.range_id_1"]
sc_3=mmcif_dict["_pdbx_struct_sheet_hbond.range_id_2"]
sc_4=mmcif_dict["_pdbx_struct_sheet_hbond.range_1_label_atom_id"]
sc_5=mmcif_dict["_pdbx_struct_sheet_hbond.range_1_label_comp_id"]
print(sc_1)
print(sc_2)
print(sc_3)
print(sc_4)
print(sc_5)
sc=mmcif_dict["_entity_poly_seq.mon_id"]
print(sc)                          #print all amino acides from an MMCIF file
y_list=mmcif_dict["_atom_site.Cartn_y"]
print(y_list)                      #get the list of the y coordinates of all atoms       
###########################################################################################
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser()
structure = parser.get_structure('2FAT', '2FAT.cif')
#child_list = parent_entity.get_list()
#parent_entity = child_entity.get_parent()
resnames = structure.get_residues()   
print (list(resnames))
######################################################################
from Bio.PDB import *
pdbl = PDBList()
pdbl.retrieve_pdb_file('4XP1')    #download file in the working directory

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('4XP1', 'xp/4xp1.cif')
struc_dict = MMCIF2Dict('xp/4xp1.cif')

ss2=struc_dict['_citation.title']
print(ss2)  #print citation title of the file
ss3=struc_dict["_struct_ref_seq.pdbx_PDB_id_code"]
print(ss3)  #print mmcif id
ss1=struc_dict["_pdbx_entity_nonpoly.name"]
print(ss1)  #print non polymer entity 
ss4=struc_dict['_chem_comp.formula']
print(ss4)  #print chemical formula of non polymer entity
################################

print(structure.header.keys())      #print headers of the file
def cleandir(obj):
    print(", ".join([a for a in dir(obj) if not a.startswith("_")]))
cleandir(structure)             #print functions that we need to use to extract info
##############################
for model in structure:
    print(f"model {model}")                           #print id of model

model = structure[0] #since we only have one model
for chain in model:
    print(f"chain {chain}, Chain ID: {chain.id}")    #print chains exist in model
    
chain_A = model['A']  
for res in chain_A:                   #print residue_name,disorder,has a certain atom or not,res_id,res_full id and is a amino acid or not that exists in chain_A
    print(f"Residue name: {res.resname}, disorder: {res.is_disordered()} ,has a certain atom or not: {res.has_id(res)}, res_id: {res.id}, res_full id: {res.get_full_id()},is a amino acid or not: {is_aa(res)}") #show all not only number id

chain_A = model['A']
for res in chain_A:                       #print residue_name and res_id that exists in chain_A
    print(f"Residue name: {res.resname}, res_id: {res.id[1]}")

chain_L = model['L']
for res in chain_L:                                            #print residue_name and res_id that exists in chain_L
    print(f"Residue name: {res.resname}, res_id: {res.id[1]}")
    
chain_H = model['H']
for res in chain_H:                                             #print residue_name and res_id that exists in chain_H
    print(f"Residue name: {res.resname}, res_id: {res.id[1]}")
    
chain_B = model['B']
for res in chain_B:                                             #print residue_name and res_id that exists in chain_B
    print(f"Residue name: {res.resname}, res_id: {res.id[1]}")
    
chain_C = model['C']
for res in chain_C:                                             #print residue_name and res_id that exists in chain_C
    print(f"Residue name: {res.resname}, res_id: {res.id[1]}")

res = chain_A[56]
print(res)         #print residue that has id=56

for atom in res:     #print atom_name,atom_fullName,atom_id,coordinate,coordinate as vector,bfactor,occupancy and disorder of the residue that has id=56
    print(f" name: {atom.name},fullName: {atom.get_fullname()} ,id: {atom.get_id()}, coordinate: {atom.get_coord()},coordinates as vector: {atom.get_vector()},bfactor: {atom.get_bfactor()},occupancy: {atom.get_occupancy()},disorder: {atom.is_disordered()}")
#########################################################################################
LDP = None
for res in structure[0].get_residues():
    if res.resname == "LDP":
        LDP = res
        break
print(LDP)                 #finding LDP residue

for res in structure[0].get_residues():
    if res.resname == "LDP":
        for atoms in res.get_atoms():
           print(atoms.name , atoms.get_coord())        #print all atoms in residue LDP

#get the distance along each axis between the two CA atoms of the two residues
res_1_CA = structure[0]['A'][56]['CA']
print(res_1_CA.coord)

res_2_CA = structure[0]['A'][327]['CA']
print(res_2_CA.coord)

diff=res_1_CA.coord - res_2_CA.coord
print(diff)

#To get positive values we square the vector and then take the square root.
import numpy as np
dist = np.sqrt(diff * diff)
print(dist)

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

# %%
