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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects if the synthesis route employs aromatic halogenation,
    specifically looking for halogenation of aromatic rings.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for aromatic halogenation reactions
            halogenation_reactions = [
                "Aromatic fluorination",
                "Aromatic chlorination",
                "Aromatic bromination",
                "Aromatic iodination",
            ]

            # First try to directly identify the reaction type
            for reaction_type in halogenation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found {reaction_type} at depth {depth}")
                    result = True
                    return

            # If specific reaction check fails, try to detect by analyzing reactants and products
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has aromatic halide
                if checker.check_fg("Aromatic halide", product):
                    # Define common halogenation reagents to exclude
                    common_reagents = [
                        "Cl2",
                        "Br2",
                        "I2",
                        "F2",
                        "NCl",
                        "NBr",
                        "NI",
                        "NF",
                        "HCl",
                        "HBr",
                        "HI",
                        "HF",
                        "ICl",
                        "IBr",
                    ]

                    # Find main reactant (excluding common reagents)
                    main_reactant = None
                    for r in reactants:
                        # Skip common reagents
                        if len(r) > 5 and not any(reagent in r for reagent in common_reagents):
                            main_reactant = r
                            break

                    if main_reactant:
                        # If main reactant doesn't have aromatic halide, it's likely a halogenation
                        if not checker.check_fg("Aromatic halide", main_reactant):
                            print(f"Found aromatic halogenation (FG analysis) at depth {depth}")
                            result = True
                            return
                        else:
                            # Count aromatic halides in both using the checker function
                            product_halide_indices = checker.get_fg_atom_indices(
                                "Aromatic halide", product
                            )
                            reactant_halide_indices = checker.get_fg_atom_indices(
                                "Aromatic halide", main_reactant
                            )

                            if len(product_halide_indices) > len(reactant_halide_indices):
                                print(
                                    f"Found aromatic halogenation (count analysis) at depth {depth}"
                                )
                                result = True
                                return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
