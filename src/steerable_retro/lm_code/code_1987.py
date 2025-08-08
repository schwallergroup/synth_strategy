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
    Detects if the synthesis includes an ester hydrolysis step.
    """
    hydrolysis_found = False

    def dfs_traverse(node):
        nonlocal hydrolysis_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")
                        if reactant_mol.HasSubstructMatch(ester_pattern):
                            # Check for carboxylic acid in product
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[#8])")
                                if product_mol.HasSubstructMatch(acid_pattern):
                                    print("Found ester hydrolysis")
                                    hydrolysis_found = True
                                    break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return hydrolysis_found
