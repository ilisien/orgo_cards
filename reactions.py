import re
from typing import List
from rdkit import Chem
from mol_handling import is_valid_smiles, iupac_to_smiles

class Reaction:
    def __init__(self, reactants: List[str], reagents: List[List[str]], products: List[str]):
        self.reactants = reactants
        self.reagents = reagents
        self.products = products

    def __repr__(self):
        return f"Reaction(reactants={self.reactants}, reagents={self.reagents}, products={self.products})"

def format_for_compound_string(compound):
    if isinstance(compound, list):
        return [format_for_compound_string(element) for element in compound]
    elif compound[0] == '"':
        return "&#8203;" + compound.replace('"',"")
    elif is_valid_smiles(compound):
        return compound
    else:
        try:
            return iupac_to_smiles(compound)
        except:
            return "error at formatting from .txt file"

def parse_reaction(line: str) -> Reaction:
    parts = line.strip().split('>')
    
    reactants = format_for_compound_string(parts[0].split(';'))
    reagents = format_for_compound_string([p.split(';') for p in parts[1].split('/')] if '/' in parts[1] else [parts[1].split(';'), []])
    products = format_for_compound_string(parts[2].split(';'))
    
    return Reaction(reactants, reagents, products)


def parse_reaction_file(filename: str) -> List[Reaction]:
    reactions = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                reactions.append(parse_reaction(line))
    return reactions


if __name__ == "__main__":
    filename = "rxn_lists/test.txt"  # Replace with actual filename
    reactions = parse_reaction_file(filename)
    for reaction in reactions:
        print(reaction)
