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
    This function detects a synthetic strategy involving the formation of a tetrahydroisoquinoline
    ring system in the middle of the synthesis.
    """
    tetrahydroisoquinoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7]=[#6][#6]1")
    amide_pattern = Chem.MolFromSmarts("[#6][#7][#6](=[#8])[#6]")

    found_thiq_formation = False

    def dfs_traverse(node):
        nonlocal found_thiq_formation

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check if product has tetrahydroisoquinoline but reactants don't
                    product_has_thiq = product.HasSubstructMatch(tetrahydroisoquinoline_pattern)
                    reactants_have_thiq = any(
                        r.HasSubstructMatch(tetrahydroisoquinoline_pattern) for r in reactants
                    )

                    # Check if reactants have amide but product doesn't (amide to imine conversion)
                    reactants_have_amide = any(
                        r.HasSubstructMatch(amide_pattern) for r in reactants
                    )

                    if product_has_thiq and not reactants_have_thiq and reactants_have_amide:
                        print(f"Found tetrahydroisoquinoline formation: {rsmi}")
                        found_thiq_formation = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_thiq_formation
