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
    This function detects a synthetic strategy involving late-stage Suzuki coupling
    between aromatic systems.
    """
    late_stage_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_suzuki

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider reactions at depth 0 or 1 (late stage)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for boronic acid in reactants
            boronic_acid_present = False
            aryl_halide_present = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                    if mol.HasSubstructMatch(boronic_pattern):
                        boronic_acid_present = True

                    aryl_halide_pattern = Chem.MolFromSmarts("[#6;a]-[#35,#53,#17]")
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True

            # Check if both boronic acid and aryl halide are present (indicative of Suzuki coupling)
            if boronic_acid_present and aryl_halide_present:
                # Verify C-C bond formation between aromatic rings
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # This is a simplified check - a more robust implementation would
                    # track specific atoms involved in the coupling
                    late_stage_suzuki = True
                    print(f"Late-stage Suzuki coupling detected at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage Suzuki coupling strategy detected: {late_stage_suzuki}")
    return late_stage_suzuki
