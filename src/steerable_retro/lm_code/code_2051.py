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
    Detects a strategy involving late-stage oxidation (alcohol to ketone).
    In retrosynthetic analysis, this appears as a reduction of ketone to alcohol.
    """
    found_late_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_oxidation

        if node["type"] == "reaction" and depth <= 2:  # Late stage (within first few reactions)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, we're looking for ketone reduction (appears as alcohol oxidation in forward direction)
                # Check if this is a reduction reaction (ketone to alcohol in retrosynthesis)
                if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                    print(f"Found ketone reduction reaction at depth {depth}")

                    # Verify reactant has ketone and product has secondary alcohol (retrosynthetic perspective)
                    reactant_has_ketone = any(
                        checker.check_fg("Ketone", reactant) for reactant in reactants
                    )
                    product_has_alcohol = checker.check_fg("Secondary alcohol", product)

                    if reactant_has_ketone and product_has_alcohol:
                        print(
                            f"Confirmed: Ketone in reactant reduced to secondary alcohol in product (retrosynthetic)"
                        )
                        found_late_oxidation = True
                        return

                # Check for oxidation reaction (will appear as reduction in retrosynthesis)
                if checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                ):
                    print(f"Found oxidation reaction at depth {depth}")

                    # In forward direction: alcohol → ketone
                    # In retrosynthesis: product has alcohol, reactant has ketone
                    reactant_has_ketone = any(
                        checker.check_fg("Ketone", reactant) for reactant in reactants
                    )
                    product_has_alcohol = checker.check_fg("Secondary alcohol", product)

                    if reactant_has_ketone and product_has_alcohol:
                        print(
                            f"Confirmed: Ketone in reactant transformed to secondary alcohol in product (retrosynthetic)"
                        )
                        found_late_oxidation = True
                        return

                # Alternative check: look for the functional group transformation directly
                # In retrosynthesis: ketone (reactant) → secondary alcohol (product)
                reactant_has_ketone = any(
                    checker.check_fg("Ketone", reactant) for reactant in reactants
                )
                product_has_alcohol = checker.check_fg("Secondary alcohol", product)

                if reactant_has_ketone and product_has_alcohol:
                    print(f"Found potential ketone to alcohol transformation at depth {depth}")
                    # Try to verify this is a reduction reaction
                    if not found_late_oxidation:  # Only check if we haven't found it yet
                        found_late_oxidation = True
                        print(f"Confirmed functional group transformation at depth {depth}")
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_oxidation
