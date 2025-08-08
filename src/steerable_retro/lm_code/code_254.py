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
    Detects if the synthesis route involves a late-stage nucleophilic substitution
    replacing a methylsulfonyl group with a thiazolidine-1,1-dioxide group.
    """
    late_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal late_substitution

        if node["type"] == "reaction" and depth <= 1:  # Very late stage
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for methylsulfonyl in reactants
                methylsulfonyl_pattern = Chem.MolFromSmarts("[#6][S](=[O])(=[O])[#6]")
                # Check for thiazolidine-1,1-dioxide in product
                thiazolidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][S]1(=[O])=[O]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(thiazolidine_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(methylsulfonyl_pattern):
                            print(
                                "Detected late-stage nucleophilic substitution with thiazolidine-1,1-dioxide"
                            )
                            late_substitution = True
                            break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_substitution
