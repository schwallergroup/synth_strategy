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
    This function detects a synthetic strategy that maintains a halogenated aromatic
    scaffold throughout the synthesis.
    """
    # Initialize tracking variables
    total_molecules = 0
    halogenated_aromatics = 0
    preserved_scaffold_reactions = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_molecules, halogenated_aromatics, preserved_scaffold_reactions, total_reactions

        if node["type"] == "mol":
            total_molecules += 1
            try:
                # Check if molecule contains aromatic halide
                if checker.check_fg("Aromatic halide", node["smiles"]):
                    halogenated_aromatics += 1
                    print(f"Found halogenated aromatic: {node['smiles']}")
            except Exception as e:
                print(f"Error processing molecule: {node['smiles']} - {str(e)}")

        elif node["type"] == "reaction":
            total_reactions += 1
            try:
                # Extract product and reactants from reaction SMILES
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product has halogenated aromatic
                    product_has_halogen = checker.check_fg("Aromatic halide", product)

                    # Check if at least one reactant has halogenated aromatic
                    reactant_has_halogen = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    )

                    # If both product and at least one reactant have halogenated aromatic,
                    # consider this as preserving the scaffold
                    if product_has_halogen and reactant_has_halogen:
                        preserved_scaffold_reactions += 1
                        print(f"Scaffold preserved in reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {str(e)}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy criteria:
    # 1. At least 50% of molecules should have halogenated aromatics
    # 2. At least 75% of reactions should preserve the halogenated aromatic scaffold
    mol_criterion = total_molecules > 0 and halogenated_aromatics / total_molecules >= 0.5
    rxn_criterion = total_reactions > 0 and preserved_scaffold_reactions / total_reactions >= 0.75

    strategy_present = mol_criterion and (total_reactions == 0 or rxn_criterion)

    print(f"Mol criterion: {mol_criterion}, Rxn criterion: {rxn_criterion or total_reactions == 0}")
    print(f"Halogenated aromatic scaffold strategy detected: {strategy_present}")
    print(f"Total molecules: {total_molecules}")
    print(f"Halogenated aromatics: {halogenated_aromatics}")
    print(f"Total reactions: {total_reactions}")
    print(f"Reactions preserving scaffold: {preserved_scaffold_reactions}")

    return strategy_present
