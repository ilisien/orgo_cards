from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchi, MolFromInchi
from rdkit.Chem import AllChem
import pubchempy as pcp

def smiles_to_iupac(smiles,synonyms_amount=3):
    """Convert a SMILES string to an IUPAC name."""
    mol = Chem.MolFromSmiles(smiles)
    inchi = MolToInchi(mol)
    compound = pcp.get_compounds(inchi, 'inchi')
    iupac_name = compound[0].iupac_name if compound else "Not found"
    compound[0].synonyms.remove(iupac_name)
    synonyms = compound[0].synonyms[:synonyms_amount] if compound else "Not found"
    return iupac_name, synonyms

def iupac_to_smiles(iupac_name):
    """Convert an IUPAC name to a SMILES string."""
    compounds = pcp.get_compounds(iupac_name, 'name')
    smiles = compounds[0].canonical_smiles if compounds else "Not found"
    return smiles

def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol)

if __name__ == "__main__":
    smiles = "CCC=O"
    print(f"SMILES: {smiles}")

    iupac_name, synonyms = smiles_to_iupac(smiles,5)
    print(f"IUPAC Name: {iupac_name}")
    print(f"Synonyms: {synonyms}")

    new_smiles = iupac_to_smiles(iupac_name)
    print(f"Reconstructed SMILES: {new_smiles}")

    draw_molecule(new_smiles).show()