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
    Detects a synthesis strategy involving cyclopropanation in the middle steps
    followed by late-stage amide formation.
    """
    # Initialize tracking variables
    has_cyclopropanation = False
    has_late_amide_formation = False

    # Track reaction depths for sequence analysis
    cyclopropanation_depth = None
    amide_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_cyclopropanation, has_late_amide_formation
        nonlocal cyclopropanation_depth, amide_formation_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for cyclopropanation
                # Method 1: Direct reaction check
                if checker.check_reaction("Diazo addition", rsmi):
                    print(f"Found cyclopropanation via Diazo addition at depth {depth}")
                    has_cyclopropanation = True
                    cyclopropanation_depth = depth
                # Method 2: Check for diazo reactant and cyclopropane product
                else:
                    diazo_present = any(
                        checker.check_fg("Diazo", reactant) for reactant in reactants
                    )
                    if diazo_present and checker.check_ring("cyclopropane", product):
                        # Verify cyclopropane wasn't already in reactants
                        reactants_have_cyclopropane = any(
                            checker.check_ring("cyclopropane", reactant) for reactant in reactants
                        )
                        if not reactants_have_cyclopropane:
                            print(f"Found cyclopropanation with diazo compound at depth {depth}")
                            has_cyclopropanation = True
                            cyclopropanation_depth = depth

                # Check for amide formation
                # Method 1: Direct reaction checks for various amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                ]

                for rxn_name in amide_formation_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Found amide formation ({rxn_name}) at depth {depth}")
                        has_late_amide_formation = True
                        amide_formation_depth = depth
                        break

                # Method 2: Check for amide product that wasn't in reactants
                if not has_late_amide_formation:
                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                    product_has_amide = any(
                        checker.check_fg(amide_type, product) for amide_type in amide_types
                    )

                    if product_has_amide:
                        # Check if reactants already had the amide
                        reactants_have_amide = any(
                            any(
                                checker.check_fg(amide_type, reactant) for amide_type in amide_types
                            )
                            for reactant in reactants
                        )

                        if not reactants_have_amide:
                            print(f"Found amide formation (FG analysis) at depth {depth}")
                            has_late_amide_formation = True
                            amide_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present:
    # 1. Both cyclopropanation and amide formation must be detected
    # 2. Cyclopropanation should occur at a greater depth than amide formation (earlier in synthesis)
    # 3. Amide formation should be late-stage (depth ≤ 1)

    strategy_present = (
        has_cyclopropanation
        and has_late_amide_formation
        and cyclopropanation_depth is not None
        and amide_formation_depth is not None
        and cyclopropanation_depth > amide_formation_depth
        and amide_formation_depth <= 1
    )

    print(f"Strategy detection result: {strategy_present}")
    print(
        f"Cyclopropanation depth: {cyclopropanation_depth}, Amide formation depth: {amide_formation_depth}"
    )

    return strategy_present
