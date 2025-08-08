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
    This function detects a strategy involving esterification using an acyl chloride.
    """
    found_acyl_chloride_esterification = False

    def dfs_traverse(node):
        nonlocal found_acyl_chloride_esterification

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Acyl chloride pattern
                    acyl_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")
                    # Phenol pattern
                    phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
                    # Ester pattern
                    ester_pattern = Chem.MolFromSmarts("[c]-[O]-[C](=[O])-[#6]")

                    # Check for acyl chloride and phenol in reactants, and ester in product
                    if (
                        reactant_mol.HasSubstructMatch(acyl_chloride_pattern)
                        and reactant_mol.HasSubstructMatch(phenol_pattern)
                        and product_mol.HasSubstructMatch(ester_pattern)
                    ):
                        found_acyl_chloride_esterification = True
                        print(f"Found acyl chloride esterification: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Acyl chloride esterification detected: {found_acyl_chloride_esterification}")
    return found_acyl_chloride_esterification
