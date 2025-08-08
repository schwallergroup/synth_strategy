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
    This function detects early-stage heterocycle formation, specifically oxazole ring formation.
    """
    heterocycle_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found

        if node["type"] == "reaction" and depth >= 2:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxazole in product but not in reactants
                oxazole_pattern = Chem.MolFromSmarts("[#6]1:[#7]:[#6]:[#8]:[#6]:1")

                has_oxazole_in_reactants = False
                for reactant in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol and r_mol.HasSubstructMatch(oxazole_pattern):
                            has_oxazole_in_reactants = True
                            break
                    except:
                        continue

                if not has_oxazole_in_reactants:
                    try:
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol and p_mol.HasSubstructMatch(oxazole_pattern):
                            print("Detected early-stage heterocycle (oxazole) formation")
                            heterocycle_formation_found = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return heterocycle_formation_found
