from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchi, MolFromInchi
from rdkit.Chem import AllChem
import pubchempy as pcp

def smiles_to_iupac(smiles):
    """Convert a SMILES string to an IUPAC name."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    inchi = MolToInchi(mol)  # Convert to InChI first
    compound = pcp.get_compounds(inchi, 'inchi')
    iupac_name = compound[0].iupac_name if compound else "Not found"
    return iupac_name

def iupac_to_smiles(iupac_name):
    """Convert an IUPAC name to a SMILES string."""
    compounds = pcp.get_compounds(iupac_name, 'name')
    smiles = compounds[0].canonical_smiles if compounds else "Not found"
    return smiles

if __name__ == "__main__":
    smiles = "CO"
    print(f"SMILES: {smiles}")

    iupac_name = smiles_to_iupac(smiles)
    print(f"IUPAC Name: {iupac_name}")

    new_smiles = iupac_to_smiles(iupac_name)
    print(f"Reconstructed SMILES: {new_smiles}")