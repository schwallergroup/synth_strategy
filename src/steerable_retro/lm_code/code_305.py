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
    Detects if the synthesis involves coupling with an aromatic fragment (like naphthalene)
    in the late stage of the synthesis.
    """
    aromatic_coupling_late = False
    late_stage_threshold = 1  # Consider depth <= 1 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_coupling_late

        if node["type"] == "reaction" and depth <= late_stage_threshold:
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for naphthalene or other complex aromatic systems in reactants
                naphthalene_pattern = Chem.MolFromSmarts("[c]1[c][c][c]2[c][c][c][c][c]2[c]1")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(naphthalene_pattern):
                        # Verify this fragment is incorporated into the product
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(naphthalene_pattern):
                            print(
                                f"Late-stage aromatic fragment coupling detected at depth {depth}"
                            )
                            aromatic_coupling_late = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return aromatic_coupling_late
