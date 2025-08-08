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
    This function detects late-stage N-hydroxyalkylation via epoxide ring opening.
    """
    has_epoxide_opening = False

    def dfs_traverse(node, depth=0):
        nonlocal has_epoxide_opening

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for epoxide in reactants
                epoxide_pattern = Chem.MolFromSmarts("[#6]1[#8][#6]1")
                alcohol_pattern = Chem.MolFromSmarts("[#8H1]")

                has_epoxide = False
                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(epoxide_pattern):
                            has_epoxide = True
                            break
                    except:
                        continue

                # Check for alcohol in product (result of epoxide opening)
                has_alcohol = False
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(alcohol_pattern):
                        has_alcohol = True
                except:
                    pass

                if has_epoxide and has_alcohol:
                    print("Found late-stage epoxide opening reaction")
                    has_epoxide_opening = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_epoxide_opening
