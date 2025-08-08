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
    Detects a linear synthesis strategy where each step involves modification of one functional group,
    with a specific sequence of functional group interconversions (e.g., nitro → amine → amide).
    """
    # Initialize tracking variables
    reactions_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Identify reaction type
                nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                amine_pattern = Chem.MolFromSmarts("[#7H2]")
                amide_pattern = Chem.MolFromSmarts("[#7H]-[#6](=[#8])")
                chloro_pattern = Chem.MolFromSmarts("[Cl]-[#6]")

                # Check for nitro reduction
                if any(
                    r.HasSubstructMatch(nitro_pattern) for r in reactants if r
                ) and product.HasSubstructMatch(amine_pattern):
                    reactions_sequence.append(("nitro_reduction", depth))
                    print(f"Found nitro reduction at depth {depth}")

                # Check for amide formation
                if any(
                    r.HasSubstructMatch(amine_pattern) for r in reactants if r
                ) and product.HasSubstructMatch(amide_pattern):
                    reactions_sequence.append(("amide_formation", depth))
                    print(f"Found amide formation at depth {depth}")

                # Check for chloride displacement
                if any(
                    r.HasSubstructMatch(chloro_pattern) for r in reactants if r
                ) and not product.HasSubstructMatch(chloro_pattern):
                    reactions_sequence.append(("chloride_displacement", depth))
                    print(f"Found chloride displacement at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (higher depth = earlier in synthesis)
    reactions_sequence.sort(key=lambda x: x[1], reverse=True)
    reaction_types = [r[0] for r in reactions_sequence]

    # Check if we have the specific sequence: chloride displacement → nitro reduction → amide formation
    expected_sequence = ["chloride_displacement", "nitro_reduction", "amide_formation"]
    sequence_match = len(reaction_types) >= 3 and all(
        exp == act for exp, act in zip(expected_sequence, reaction_types[:3])
    )

    print(f"Reaction sequence: {reaction_types}")
    print(f"Linear functional group interconversion strategy detected: {sequence_match}")
    return sequence_match
