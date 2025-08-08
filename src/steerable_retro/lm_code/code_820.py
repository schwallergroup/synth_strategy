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
    This function detects if the synthetic route involves a late-stage acylation
    with an acryloyl chloride as the final step.
    """
    final_step_is_acylation = False

    def dfs_traverse(node):
        nonlocal final_step_is_acylation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Check if this is the root reaction (depth 0)
            if len(node.get("children", [])) == 2:  # Typical reaction has 2 children
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if one reactant is an acyl chloride (specifically acryloyl chloride)
                acyl_chloride_pattern = Chem.MolFromSmarts("C=CC(=O)Cl")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(acyl_chloride_pattern):
                            # Check if product has an amide bond
                            product_mol = Chem.MolFromSmiles(product)
                            amide_pattern = Chem.MolFromSmarts("C=CC(=O)N")
                            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                                final_step_is_acylation = True
                                print("Detected late-stage acylation with acryloyl chloride")
                    except:
                        continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return final_step_is_acylation
