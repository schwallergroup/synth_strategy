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
    Detects if the synthesis involves nitration in a late stage (low depth).
    """
    late_nitration_found = False

    def dfs_traverse(node):
        nonlocal late_nitration_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction: {rsmi}")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth information with a default value
            depth_match = re.search(r"Depth: (\d+)", node.get("metadata", {}).get("ID", ""))
            depth = int(depth_match.group(1)) if depth_match else 0

            # Only consider late stage reactions (depth 0, 1, or 2)
            if depth <= 2:
                print(f"Checking late-stage reaction at depth {depth}")

                # Check for nitration reactions using the checker function
                nitration_reaction_types = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Non-aromatic nitration with HNO3",
                ]

                is_nitration = False
                for rxn_type in nitration_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_nitration = True
                        print(f"Found nitration reaction type: {rxn_type}")
                        break

                # If not identified as a standard nitration reaction, check for nitro group appearance
                if not is_nitration:
                    # Check if nitro group appears in product but not in reactants
                    has_nitro_in_product = checker.check_fg("Nitro group", product)
                    has_nitro_in_reactants = any(
                        checker.check_fg("Nitro group", reactant) for reactant in reactants
                    )

                    if has_nitro_in_product and not has_nitro_in_reactants:
                        # Additional check to confirm this is actually a nitration reaction
                        # and not just a reaction where a nitro group appears for other reasons
                        is_nitration = True
                        print(f"Found nitration by nitro group appearance at depth {depth}")

                if is_nitration:
                    late_nitration_found = True
                    print(f"Found late-stage nitration at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if not late_nitration_found:
        print("No late-stage nitration found in the route")

    return late_nitration_found
