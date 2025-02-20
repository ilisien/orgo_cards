from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.inchi import MolToInchi, MolFromInchi
from rdkit.Chem import AllChem

def smiles_to_iupac(smiles):
    """Convert a SMILES string to an IUPAC name."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    inchi = MolToInchi(mol)  # Convert to InChI first
    mol_from_inchi = MolFromInchi(inchi)  # Convert back to molecule
    iupac_name = Chem.MolToIUPACName(mol_from_inchi)  # Get IUPAC name
    return iupac_name

def iupac_to_smiles(iupac_name):
    """Convert an IUPAC name to a SMILES string."""
    mol = Chem.MolFromIUPACName(iupac_name)
    if mol is None:
        raise ValueError("Invalid IUPAC name.")
    smiles = Chem.MolToSmiles(mol)
    return smiles

def print_ascii_representation(smiles):
    """Print an ASCII representation of a molecule from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    print("ASCII Representation:")
    print(Chem.MolToMolBlock(mol))  # Print MolBlock (includes ASCII representation)

# Example usage
smiles = "C=O"  # Formaldehyde
print(f"SMILES: {smiles}")

# Convert SMILES to IUPAC
iupac_name = smiles_to_iupac(smiles)
print(f"IUPAC Name: {iupac_name}")

# Convert IUPAC back to SMILES
new_smiles = iupac_to_smiles(iupac_name)
print(f"Reconstructed SMILES: {new_smiles}")

# Print ASCII representation
print_ascii_representation(smiles)