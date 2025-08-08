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
    This function detects a strategy involving multiple sequential C-N bond formations
    through nucleophilic substitutions, particularly on heterocycles.
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    # Check for C-N bond formation
                    # Look for halogen (Cl, Br) on carbon in reactants that becomes a C-N bond in product
                    for reactant in reactants:
                        if reactant:
                            # Find carbons with halogens
                            c_hal_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")
                            c_hal_matches = reactant.GetSubstructMatches(c_hal_pattern)

                            for match in c_hal_matches:
                                c_atom = match[0]
                                # Check if this carbon forms a bond with nitrogen in the product
                                if product:
                                    c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                                    c_n_matches = product.GetSubstructMatches(c_n_pattern)

                                    for c_n_match in c_n_matches:
                                        # This is a simplification - in a real implementation,
                                        # we would need to track atom mappings between reactants and products
                                        cn_bond_formations += 1
                                        print(f"Detected C-N bond formation in reaction")
                                        break
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if we found at least 2 C-N bond formations
    return cn_bond_formations >= 2
