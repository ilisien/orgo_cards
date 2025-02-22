import re, os, time
from typing import List
from rdkit import Chem
from mol_handling import is_valid_smiles, iupac_to_smiles, draw_molecule, smiles_to_iupac
from flask import url_for

class Reaction:
    def __init__(self, reactants: List[str], reagents: List[List[str]], products: List[str]):
        self.reactants = reactants
        self.reagents = reagents
        self.products = products
        self.full_rxn_list = [reactants,reagents,products]

    def __repr__(self):
        return f"Reaction(reactants={self.reactants}, reagents={self.reagents}, products={self.products})"
    
    def load_compound_images(self,rxn=None):
        if rxn is None:
            self.load_compound_images(self.full_rxn_list)
        elif isinstance(rxn,list):
            for r in rxn:
                self.load_compound_images(r)
        elif is_valid_smiles(rxn):
            iupac = smiles_to_iupac(rxn)[0].replace(" ","_")
            if not os.path.exists(f"static/images/{iupac}.png"):
                draw_molecule(rxn).save(f"static/images/{iupac}.png")
            time.sleep(1)
    
    def get_html_components(self,hide):
        return {
            "hide" : hide,
            "image_reactants" : [smiles_to_iupac(reactant)[0] for reactant in self.reactants if is_valid_smiles(reactant)],
            "nonimage_reactants" : [reactant for reactant in self.reactants if not is_valid_smiles(reactant)],
            "top_reagents" : self.reagents[0],
            "bottom_reagents" : self.reagents[1],
            "image_products" : [smiles_to_iupac(product)[0] for product in self.products if is_valid_smiles(product)],
            "nonimage_products" : [product for product in self.products if not is_valid_smiles(product)],
        }

    def all_html_components(self):
        self.load_compound_images()
        return [self.get_html_components(section) for section in [0,1,2]]

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
    html_components = []
    for reaction in reactions:
        html_components.extend(reaction.all_html_components())
    print(html_components)