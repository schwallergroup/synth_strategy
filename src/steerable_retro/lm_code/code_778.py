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
    This function detects if the final step (or one of the last steps)
    in the synthesis is an amide formation.
    """
    # Track depth of amide formation
    amide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]

                    if product_mol and reactant_mols:
                        # Check for amide formation
                        amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                        amine_pattern = Chem.MolFromSmarts("[N;!$(NC=O)]")

                        # Check if product has amide
                        if product_mol.HasSubstructMatch(amide_pattern):
                            # Check if reactants have carboxylic acid and amine
                            has_acid = any(
                                mol.HasSubstructMatch(carboxylic_acid_pattern)
                                for mol in reactant_mols
                            )
                            has_amine = any(
                                mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                            )

                            if has_acid and has_amine:
                                amide_formation_depth = depth
                except:
                    print("Error processing molecule in late_stage_amide_formation")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if amide formation is in the last 30% of steps
    if amide_formation_depth is not None and amide_formation_depth <= max_depth * 0.3:
        print(
            f"Detected late-stage amide formation at depth {amide_formation_depth} of {max_depth}"
        )
        return True
    return False
