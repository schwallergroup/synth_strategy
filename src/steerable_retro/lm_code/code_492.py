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
    Detects if the synthetic route uses click chemistry (azide-alkyne cycloaddition)
    to form a triazole ring in the final step of synthesis.
    """
    # Track if we found click chemistry in a late stage
    found_click_chemistry = False

    # Get the target molecule SMILES (root of the tree)
    target_molecule_smiles = route["smiles"]
    print(f"Target molecule SMILES: {target_molecule_smiles}")

    def dfs_traverse(node, depth=0):
        nonlocal found_click_chemistry

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a late-stage reaction (depth <= 1)
                is_late_stage = depth <= 1
                print(f"Is late-stage reaction: {is_late_stage}")

                # Check if the product is the target molecule or one step away
                is_target_or_close = (product_smiles == target_molecule_smiles) or (depth <= 1)
                print(f"Is target molecule or close: {is_target_or_close}")

                if is_late_stage:
                    # Check for azide in reactants
                    has_azide = any(checker.check_fg("Azide", r) for r in reactants_smiles)
                    print(f"Has azide in reactants: {has_azide}")

                    # Check for alkyne in reactants
                    has_alkyne = any(checker.check_fg("Alkyne", r) for r in reactants_smiles)
                    print(f"Has alkyne in reactants: {has_alkyne}")

                    # Check for triazole in product
                    has_triazole_in_product = checker.check_ring("triazole", product_smiles)
                    print(f"Has triazole in product: {has_triazole_in_product}")

                    # Check if any reactant already has triazole
                    has_triazole_in_reactants = any(
                        checker.check_ring("triazole", r) for r in reactants_smiles
                    )
                    print(f"Has triazole in reactants: {has_triazole_in_reactants}")

                    # Check if this is a Huisgen cycloaddition reaction
                    is_huisgen_reaction = (
                        checker.check_reaction(
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction("Huisgen_Cu-catalyzed_1,4-subst", rsmi)
                        or checker.check_reaction("Huisgen_Ru-catalyzed_1,5_subst", rsmi)
                        or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                    )
                    print(f"Is Huisgen cycloaddition: {is_huisgen_reaction}")

                    # If we have azide and alkyne in reactants, and triazole in product but not in reactants,
                    # or if it's explicitly a Huisgen cycloaddition reaction, it's click chemistry
                    if (
                        has_azide
                        and has_alkyne
                        and has_triazole_in_product
                        and not has_triazole_in_reactants
                    ) or is_huisgen_reaction:
                        found_click_chemistry = True
                        print(f"Found click chemistry triazole formation at depth {depth}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {found_click_chemistry}")

    return found_click_chemistry
