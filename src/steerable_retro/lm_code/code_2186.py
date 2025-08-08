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
    This function detects if acylation with oxalyl chloride is used in the synthesis.
    """
    oxalyl_chloride_acylation_detected = False

    def dfs_traverse(node):
        nonlocal oxalyl_chloride_acylation_detected
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for oxalyl chloride in reactants
                oxalyl_chloride_pattern = Chem.MolFromSmarts("ClC(=O)C(=O)Cl")

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(oxalyl_chloride_pattern):
                        oxalyl_chloride_acylation_detected = True
                        print(f"Acylation with oxalyl chloride detected in reaction: {rsmi}")
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Acylation with oxalyl chloride detected: {oxalyl_chloride_acylation_detected}")
    return oxalyl_chloride_acylation_detected
