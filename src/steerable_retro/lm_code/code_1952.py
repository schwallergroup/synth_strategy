#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis involves attaching a small fragment (like bromoacetic acid ester)
    to a larger fragment.
    """
    small_fragment_attachment_found = False

    def dfs_traverse(node, depth=0):
        nonlocal small_fragment_attachment_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            if len(reactants_smiles) >= 2:
                # Check sizes of reactants
                reactant_sizes = []
                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            reactant_sizes.append(mol.GetNumAtoms())
                    except:
                        continue

                # Check if we have at least one small fragment (< 10 atoms) and one larger fragment (> 15 atoms)
                if reactant_sizes and min(reactant_sizes) < 10 and max(reactant_sizes) > 15:
                    # Check if product is larger than the largest reactant
                    try:
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol and product_mol.GetNumAtoms() > max(reactant_sizes):
                            small_fragment_attachment_found = True
                            print(f"Small fragment attachment detected at depth {depth}")
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return small_fragment_attachment_found
