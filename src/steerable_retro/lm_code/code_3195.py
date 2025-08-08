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
    Detects a linear synthetic strategy starting from a Weinreb amide,
    proceeding through a ketone and α-bromoketone, to form a thiazole.
    """
    # Track the sequence of transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and all(reactants):
                    # Define patterns
                    weinreb_pattern = Chem.MolFromSmarts("C(=O)N(C)OC")
                    ketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6]")
                    bromoketone_pattern = Chem.MolFromSmarts("[#6]-C(=O)-C-[Br]")
                    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                    # Check for specific transformations
                    if any(
                        r.HasSubstructMatch(weinreb_pattern) for r in reactants
                    ) and product.HasSubstructMatch(ketone_pattern):
                        transformations.append(("weinreb_to_ketone", depth))
                        print(f"Detected Weinreb amide to ketone at depth {depth}")

                    if any(
                        r.HasSubstructMatch(ketone_pattern) for r in reactants
                    ) and product.HasSubstructMatch(bromoketone_pattern):
                        transformations.append(("ketone_to_bromoketone", depth))
                        print(f"Detected ketone to α-bromoketone at depth {depth}")

                    if any(
                        r.HasSubstructMatch(bromoketone_pattern) for r in reactants
                    ) and product.HasSubstructMatch(thiazole_pattern):
                        transformations.append(("bromoketone_to_thiazole", depth))
                        print(f"Detected α-bromoketone to thiazole at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all three transformations are present and in the correct order
    transformation_types = [t[0] for t in transformations]

    has_weinreb_to_ketone = "weinreb_to_ketone" in transformation_types
    has_ketone_to_bromoketone = "ketone_to_bromoketone" in transformation_types
    has_bromoketone_to_thiazole = "bromoketone_to_thiazole" in transformation_types

    # Check if the transformations are in the correct order by depth
    correct_order = True
    if has_weinreb_to_ketone and has_ketone_to_bromoketone and has_bromoketone_to_thiazole:
        depths = {t[0]: t[1] for t in transformations}
        if not (
            depths["weinreb_to_ketone"]
            > depths["ketone_to_bromoketone"]
            > depths["bromoketone_to_thiazole"]
        ):
            correct_order = False

    strategy_detected = (
        has_weinreb_to_ketone
        and has_ketone_to_bromoketone
        and has_bromoketone_to_thiazole
        and correct_order
    )

    if strategy_detected:
        print("Complete Weinreb to thiazole linear strategy detected in correct order")

    return strategy_detected
