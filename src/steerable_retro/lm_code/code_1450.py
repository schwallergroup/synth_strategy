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
    This function detects a strategy involving amide formation in the final step
    of the synthesis (depth 0).
    """
    late_stage_amide_formation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_amide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Try to extract depth from various sources
            depth = None
            if "Depth: " in rsmi:
                try:
                    depth = int(rsmi.split("Depth: ")[1].split()[0])
                except (ValueError, IndexError):
                    pass

            # If depth not found in rsmi, use current_depth from traversal
            if depth is None:
                depth = current_depth

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for amide coupling reactions directly
                    is_amide_coupling = (
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    )

                    # If direct reaction check fails, check for functional group changes
                    if not is_amide_coupling:
                        # Check reactants for amine and carboxylic acid
                        reactants_have_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )

                        reactants_have_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )

                        # Check product for amide
                        product_has_amide = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        )

                        is_amide_coupling = (
                            reactants_have_amine and reactants_have_acid and product_has_amide
                        )

                    if is_amide_coupling:
                        print(f"Late-stage amide formation detected at depth {depth}")
                        late_stage_amide_formation = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_amide_formation
