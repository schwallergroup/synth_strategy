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
    This function detects a synthetic strategy involving late-stage dehalogenation
    (dehalogenation in the final step of the synthesis).
    """
    # Track if dehalogenation occurs at the final step
    has_late_dehalogenation = False
    first_reaction_found = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_late_dehalogenation, first_reaction_found

        # Add depth attribute to the node
        node["depth"] = current_depth

        if node["type"] == "reaction":
            # Check if this is the first reaction node we encounter (final synthetic step)
            if not first_reaction_found:
                first_reaction_found = True

                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Checking final step reaction: {rsmi}")

                    # Check if this is a dehalogenation reaction using reaction checkers
                    if checker.check_reaction("Dehalogenation", rsmi) or checker.check_reaction(
                        "Aromatic dehalogenation", rsmi
                    ):
                        print(f"Detected dehalogenation reaction in final step: {rsmi}")
                        has_late_dehalogenation = True
                    else:
                        # As a fallback, check for halogen reduction manually
                        reactants_smiles = rsmi.split(">")[0].split(".")
                        product_smiles = rsmi.split(">")[-1]

                        try:
                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if not product_mol:
                                return

                            # Count halogens in product
                            product_halogen_count = sum(
                                1
                                for atom in product_mol.GetAtoms()
                                if atom.GetSymbol() in ["F", "Cl", "Br", "I"]
                            )

                            # Check each reactant for halogens
                            for reactant_smiles in reactants_smiles:
                                if not reactant_smiles:
                                    continue

                                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                                if not reactant_mol:
                                    continue

                                # Count halogens in reactant
                                reactant_halogen_count = sum(
                                    1
                                    for atom in reactant_mol.GetAtoms()
                                    if atom.GetSymbol() in ["F", "Cl", "Br", "I"]
                                )

                                # Also check for halogen functional groups to be thorough
                                has_halogen_fg = False
                                for halogen_fg in [
                                    "Primary halide",
                                    "Secondary halide",
                                    "Tertiary halide",
                                    "Aromatic halide",
                                    "Alkenyl halide",
                                    "Haloalkyne",
                                ]:
                                    if checker.check_fg(halogen_fg, reactant_smiles):
                                        has_halogen_fg = True
                                        print(f"Found {halogen_fg} in reactant: {reactant_smiles}")
                                        break

                                # If reactant has more halogens than product, it's a dehalogenation
                                if (reactant_halogen_count > product_halogen_count) or (
                                    has_halogen_fg and reactant_halogen_count > 0
                                ):
                                    print(f"Detected halogen reduction in final step: {rsmi}")
                                    print(
                                        f"Reactant halogen count: {reactant_halogen_count}, Product halogen count: {product_halogen_count}"
                                    )
                                    has_late_dehalogenation = True
                                    break

                        except Exception as e:
                            print(f"Error processing reaction: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage dehalogenation strategy detected: {has_late_dehalogenation}")
    return has_late_dehalogenation
