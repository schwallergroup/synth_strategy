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
    This function detects a strategy involving early-stage C-H functionalization
    to introduce a reactive functional group (like aldehyde).
    """
    # Initialize tracking variable
    has_early_ch_functionalization = False

    # List of functional groups commonly introduced by C-H functionalization
    target_functional_groups = [
        "Aldehyde",
        "Ketone",
        "Carboxylic acid",
        "Primary halide",
        "Secondary halide",
        "Aromatic halide",
        "Boronic acid",
        "Boronic ester",
    ]

    # C-H functionalization reaction types
    ch_functionalization_reactions = [
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "Friedel-Crafts acylation",
        "Directed ortho metalation of arenes",
        "Minisci (para)",
        "Minisci (ortho)",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_early_ch_functionalization

        if node["type"] == "reaction" and depth >= 2:  # Early stage (high depth)
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a C-H functionalization reaction
                is_ch_functionalization = False
                for reaction_type in ch_functionalization_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_ch_functionalization = True
                        print(f"Detected {reaction_type} reaction at depth {depth}")
                        break

                # If not directly identified, check for functional group introduction
                if not is_ch_functionalization:
                    for fg in target_functional_groups:
                        # Check if functional group is in product but not in any reactant
                        fg_in_product = checker.check_fg(fg, product)
                        fg_in_reactants = any(checker.check_fg(fg, r) for r in reactants)

                        if fg_in_product and not fg_in_reactants:
                            print(f"Detected introduction of {fg} at depth {depth}")
                            # Additional check for aromatic C-H functionalization
                            product_mol = Chem.MolFromSmiles(product)
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and product_mol:
                                    # Check if the reaction involves an aromatic system
                                    if any(
                                        atom.GetIsAromatic() for atom in reactant_mol.GetAtoms()
                                    ):
                                        is_ch_functionalization = True
                                        print(
                                            f"Detected aromatic C-H functionalization introducing {fg} at depth {depth}"
                                        )
                                        break

                if is_ch_functionalization:
                    has_early_ch_functionalization = True
                    print(f"Confirmed early C-H functionalization at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_early_ch_functionalization:
        print("Detected early C-H functionalization strategy")
    else:
        print("Early C-H functionalization strategy not detected")

    return has_early_ch_functionalization
