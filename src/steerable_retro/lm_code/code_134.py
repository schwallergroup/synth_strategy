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
    This function detects a synthetic strategy involving thione (C=S) to carbonyl (C=O) conversion
    as a key intermediate step.
    """
    thione_to_carbonyl_detected = False

    def dfs_traverse(node):
        nonlocal thione_to_carbonyl_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant contains thione and product contains carbonyl
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                thione_pattern = Chem.MolFromSmarts("[#6]=[#16]")
                carbonyl_pattern = Chem.MolFromSmarts("[#6]=[#8]")

                reactant_has_thione = any(
                    mol is not None and mol.HasSubstructMatch(thione_pattern)
                    for mol in reactant_mols
                )
                product_has_carbonyl = product_mol is not None and product_mol.HasSubstructMatch(
                    carbonyl_pattern
                )

                if reactant_has_thione and product_has_carbonyl:
                    print("Found thione to carbonyl conversion")
                    thione_to_carbonyl_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thione_to_carbonyl_detected
