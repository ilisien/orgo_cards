from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchi, MolFromInchi
from rdkit.Chem import AllChem
import pubchempy as pcp
from PIL import Image
import os, time

def is_valid_smiles(smiles_string):
    # unfortunately rdkit doesn't return an actual error here, so a try/except block doesn't work
    RDLogger.DisableLog('rdApp.*')
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return False
    else:
        return True

def smiles_to_iupac(smiles,synonyms_amount=0,):
    """Convert a SMILES string to an IUPAC name."""
    fast = {}
    if os.path.exists("smiles_to_iupac.txt"):
        with open("smiles_to_iupac.txt",'r') as file:
            for line in file:
                key, value = line.split(";")
                fast[key.strip()] = value.strip()
        if smiles in fast.keys():
            print(fast[smiles])
            return fast[smiles], None
    mol = Chem.MolFromSmiles(smiles)
    inchi = MolToInchi(mol)
    compound = pcp.get_compounds(inchi, 'inchi')
    iupac_name = compound[0].iupac_name if compound else "Not found"
    if iupac_name in compound[0].synonyms:
        compound[0].synonyms.remove(iupac_name)
    synonyms = compound[0].synonyms[:synonyms_amount] if compound else "Not found"
    if os.path.exists("smiles_to_iupac.txt") and smiles not in fast.keys():
        with open("smiles_to_iupac.txt","a") as f:
            print(f"{smiles};{iupac_name}", file=f)
    time.sleep(0.2) # don't overwhelm pubchem
    return iupac_name, synonyms

def iupac_to_smiles(iupac_name):
    """Convert an IUPAC name to a SMILES string."""
    compounds = pcp.get_compounds(iupac_name, 'name')
    smiles = compounds[0].canonical_smiles if compounds else "Not found"
    return smiles

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    image =  Draw.MolToImage(mol,useBWAtomPalette=True)
    bbox = image.convert("L").point(lambda x: 0 if x > 250 else 255).getbbox()
    cropped_image = image.crop(bbox)
    return cropped_image

if __name__ == "__main__":
    #smiles = "CCCO"
    #print(f"SMILES: {smiles}")

    #iupac_name, synonyms = smiles_to_iupac(smiles,5)
    #print(f"IUPAC Name: {iupac_name}")
    #print(f"Synonyms: {synonyms}")

    #new_smiles = iupac_to_smiles(iupac_name)
    #print(f"Reconstructed SMILES: {new_smiles}")

    print(iupac_to_smiles("Br2"))
    draw_molecule("CC1(C)CC1[X]").save("test.png")