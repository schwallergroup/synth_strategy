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
    Detects if the synthetic route follows a linear strategy with
    a late-stage ring formation (specifically triazole formation).
    """
    # Track reaction types and reactant counts at each depth
    reaction_types = {}
    reactant_counts = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Count number of significant reactants (excluding reagents)
            significant_reactants = 0
            for r in reactants_smiles:
                # Use atom count as better heuristic for significant molecules
                mol = Chem.MolFromSmiles(r)
                if mol and mol.GetNumAtoms() > 3:  # More than 3 atoms is significant
                    significant_reactants += 1

            reactant_counts[depth] = significant_reactants
            print(f"Significant reactants at depth {depth}: {significant_reactants}")

            # Initialize reaction type as unknown
            reaction_type = "unknown"

            # Check for triazole formation using checker function
            has_triazole_in_product = checker.check_ring("triazole", product_smiles)

            has_triazole_in_reactants = any(
                checker.check_ring("triazole", r) for r in reactants_smiles
            )

            # Check for azide and alkyne in reactants
            has_azide_in_reactants = any(checker.check_fg("Azide", r) for r in reactants_smiles)
            has_alkyne_in_reactants = any(checker.check_fg("Alkyne", r) for r in reactants_smiles)

            # Check if this is a click chemistry reaction (any triazole-forming reaction)
            is_click_reaction = (
                checker.check_reaction("Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen_Cu-catalyzed_1,4-subst", rsmi)
                or checker.check_reaction("Huisgen_Ru-catalyzed_1,5_subst", rsmi)
                or checker.check_reaction("Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi)
                or checker.check_reaction("Huisgen_disubst-alkyne", rsmi)
            )

            print(
                f"Triazole in product: {has_triazole_in_product}, Triazole in reactants: {has_triazole_in_reactants}"
            )
            print(
                f"Azide in reactants: {has_azide_in_reactants}, Alkyne in reactants: {has_alkyne_in_reactants}"
            )
            print(f"Is click reaction: {is_click_reaction}")

            # Identify triazole formation based on structural changes and reactants
            if has_triazole_in_product and not has_triazole_in_reactants:
                if is_click_reaction or (has_azide_in_reactants and has_alkyne_in_reactants):
                    reaction_type = "triazole_formation"
                    print(f"Found triazole formation at depth {depth}")
                # Even if we don't detect specific reactants or reaction type,
                # the appearance of a triazole is strong evidence
                elif has_azide_in_reactants or has_alkyne_in_reactants:
                    reaction_type = "triazole_formation"
                    print(
                        f"Found likely triazole formation at depth {depth} (detected azide or alkyne)"
                    )

            # Store reaction type at this depth
            reaction_types[depth] = reaction_type

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear synthesis
    is_linear = True
    for depth, count in reactant_counts.items():
        if count > 2:  # More than 2 significant reactants suggests convergent synthesis
            is_linear = False
            print(f"Non-linear step found at depth {depth} with {count} reactants")
            break

    # Check for late-stage triazole formation (depth 0, 1, or 2)
    max_late_stage_depth = min(2, max(reaction_types.keys()) if reaction_types else 0)
    has_late_stage_triazole = any(
        reaction_types.get(d) == "triazole_formation" for d in range(max_late_stage_depth + 1)
    )

    print(f"Is linear synthesis: {is_linear}")
    print(f"Has late-stage triazole formation: {has_late_stage_triazole}")
    print(f"Reaction types by depth: {reaction_types}")

    if is_linear and has_late_stage_triazole:
        print("Found linear synthesis with late-stage triazole formation")
        return True
    return False
