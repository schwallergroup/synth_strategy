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
    Detects if the synthesis involves an α-bromoketone intermediate that is later used
    for heterocycle formation.
    """
    has_bromoketone_intermediate = False
    bromoketone_used_for_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal has_bromoketone_intermediate, bromoketone_used_for_heterocycle

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and all(reactants):
                    bromoketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-C-[Br]")
                    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                    # Check if this reaction produces an α-bromoketone
                    if product.HasSubstructMatch(bromoketone_pattern) and not any(
                        r.HasSubstructMatch(bromoketone_pattern) for r in reactants
                    ):
                        print(f"Detected α-bromoketone formation at depth {depth}")
                        has_bromoketone_intermediate = True

                    # Check if this reaction uses an α-bromoketone to form a heterocycle
                    if any(
                        r.HasSubstructMatch(bromoketone_pattern) for r in reactants
                    ) and product.HasSubstructMatch(thiazole_pattern):
                        print(
                            f"Detected α-bromoketone used for heterocycle formation at depth {depth}"
                        )
                        bromoketone_used_for_heterocycle = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is detected if both conditions are met
    strategy_detected = has_bromoketone_intermediate and bromoketone_used_for_heterocycle
    if strategy_detected:
        print("Complete α-bromoketone intermediate strategy detected")

    return strategy_detected
