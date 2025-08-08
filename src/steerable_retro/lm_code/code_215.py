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
    This function detects if the synthetic route involves a late-stage sulfonamide formation
    (N-S bond formation in the final step).
    """
    final_reaction_has_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal final_reaction_has_sulfonamide

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Root node should be a molecule
        if depth == 0 and node["type"] != "mol":
            print("Root node is not a molecule, unexpected route structure")
            return

        # Check the first reaction (late-stage)
        if node["type"] == "reaction" and depth == 1:  # First reaction (depth 1)
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction SMILES at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if the reaction is a sulfonamide synthesis
                sulfonamide_reactions = [
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                ]

                for rxn_type in sulfonamide_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected late-stage sulfonamide formation via {rxn_type}")
                        final_reaction_has_sulfonamide = True
                        return

                # If specific reaction check fails, check for functional groups
                # Check if one reactant contains sulfonyl chloride
                has_sulfonyl_chloride = False
                for reactant in reactants:
                    if checker.check_fg("Sulfonyl halide", reactant):
                        print(f"Found sulfonyl halide in reactant: {reactant}")
                        has_sulfonyl_chloride = True
                        break

                # Check if product contains sulfonamide
                has_sulfonamide_product = checker.check_fg("Sulfonamide", product)
                if has_sulfonamide_product:
                    print(f"Found sulfonamide in product: {product}")

                # Check if any reactant already contains sulfonamide
                has_sulfonamide_reactant = False
                for reactant in reactants:
                    if checker.check_fg("Sulfonamide", reactant):
                        print(f"Found sulfonamide in reactant: {reactant}")
                        has_sulfonamide_reactant = True
                        break

                # Verify sulfonamide formation: sulfonyl chloride in reactants,
                # sulfonamide in product, and not in reactants
                if (
                    has_sulfonyl_chloride
                    and has_sulfonamide_product
                    and not has_sulfonamide_reactant
                ):
                    print(
                        f"Detected late-stage sulfonamide formation via functional group analysis"
                    )
                    final_reaction_has_sulfonamide = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {final_reaction_has_sulfonamide}")
    return final_reaction_has_sulfonamide
