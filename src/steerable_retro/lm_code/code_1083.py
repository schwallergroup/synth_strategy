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
    This function detects if the synthesis follows a linear strategy
    (each step builds on a single precursor without convergent steps).

    A linear synthesis is characterized by:
    1. Each reaction node has exactly one non-in_stock molecule child
    2. The synthesis tree forms a single path from target to starting materials
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Skip further checks if we already determined it's not linear
        if not is_linear:
            return

        indent = "  " * depth

        if node["type"] == "reaction":
            print(f"{indent}Checking reaction node at depth {depth}")

            # Count non-in_stock molecule children
            non_stock_mol_children = 0
            non_stock_smiles = []
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_mol_children += 1
                    non_stock_smiles.append(child["smiles"])

            # In a linear synthesis, each reaction should have at most one
            # non-in_stock molecule child (the main intermediate)
            if non_stock_mol_children > 1:
                print(
                    f"{indent}Found {non_stock_mol_children} non-stock molecule children - not linear"
                )
                is_linear = False
                return

            # Also check if the reaction itself is convergent by examining reactants
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for inherently convergent reaction types
                if (
                    checker.check_reaction("Suzuki", rsmi)
                    or checker.check_reaction("Diels-Alder", rsmi)
                    or checker.check_reaction("Wittig", rsmi)
                    or checker.check_reaction("Heck", rsmi)
                    or checker.check_reaction("Sonogashira", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                ):
                    print(f"{indent}Detected inherently convergent reaction type")
                    is_linear = False
                    return

                # Count significant reactants (excluding small reagents)
                significant_reactants = []
                for reactant in reactants:
                    if reactant and Chem.MolFromSmiles(reactant):
                        mol = Chem.MolFromSmiles(reactant)
                        # Consider reactants with more than 5 heavy atoms as significant
                        if mol and mol.GetNumHeavyAtoms() > 5:
                            significant_reactants.append(reactant)

                if len(significant_reactants) > 1:
                    print(
                        f"{indent}Detected convergent reaction with {len(significant_reactants)} significant reactants"
                    )
                    for r in significant_reactants:
                        print(f"{indent}  - {r}")

                    # Check if all significant reactants are in-stock
                    # If they are, this doesn't violate linearity from a synthesis planning perspective
                    all_significant_in_stock = True
                    for child in node.get("children", []):
                        if (
                            child["type"] == "mol"
                            and child["smiles"] in significant_reactants
                            and not child.get("in_stock", False)
                        ):
                            all_significant_in_stock = False
                            break

                    if not all_significant_in_stock:
                        is_linear = False
                        return

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if is_linear:
        print("Synthesis follows a linear strategy")
    else:
        print("Synthesis follows a convergent strategy")

    return is_linear
