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
    This function detects if the route includes lactam formation via nitro reduction.

    The process involves:
    1. Reduction of a nitro group to an amine
    2. Intramolecular cyclization to form a lactam ring
    """
    # Track both nitro reduction and subsequent lactam formation
    nitro_reduction_nodes = []
    lactam_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal lactam_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                # Extract reactants and product using the correct format
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for nitro reduction
                has_nitro_reactant = checker.check_fg("Nitro group", reactants_part)
                has_amine_product = checker.check_fg(
                    "Primary amine", product_part
                ) or checker.check_fg("Secondary amine", product_part)
                nitro_reduced = has_nitro_reactant and not checker.check_fg(
                    "Nitro group", product_part
                )
                is_reduction_reaction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                if (
                    has_nitro_reactant
                    and has_amine_product
                    and nitro_reduced
                    and is_reduction_reaction
                ):
                    print(f"Detected nitro reduction at depth {depth}: {rsmi}")
                    nitro_reduction_nodes.append((node, depth))

                # Check for lactam formation
                has_amide_product = (
                    checker.check_fg("Primary amide", product_part)
                    or checker.check_fg("Secondary amide", product_part)
                    or checker.check_fg("Tertiary amide", product_part)
                )

                # Check if product contains a ring
                product_mol = Chem.MolFromSmiles(product_part)
                has_ring_in_product = False
                if product_mol:
                    ring_info = product_mol.GetRingInfo()
                    has_ring_in_product = ring_info.NumRings() > 0

                # Check for direct lactam formation via nitro reduction in a single step
                if (
                    has_nitro_reactant
                    and has_amide_product
                    and has_ring_in_product
                    and nitro_reduced
                ):
                    print(
                        f"Detected direct lactam formation via nitro reduction at depth {depth}: {rsmi}"
                    )
                    lactam_formation_found = True

                # Check for lactam formation in a cyclization step after nitro reduction
                elif has_amide_product and has_ring_in_product and not has_nitro_reactant:
                    # Check if this step follows a nitro reduction step
                    for nitro_node, nitro_depth in nitro_reduction_nodes:
                        # If this step is later in the synthesis (lower depth) than a nitro reduction
                        if depth < nitro_depth:
                            print(
                                f"Detected lactam formation following nitro reduction at depth {depth}: {rsmi}"
                            )
                            print(f"  Related to nitro reduction at depth {nitro_depth}")
                            lactam_formation_found = True
                            break

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return lactam_formation_found
