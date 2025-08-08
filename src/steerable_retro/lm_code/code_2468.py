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
    This function detects a synthetic strategy involving epoxide formation,
    aminoethanol coupling, and morpholine ring formation/opening.
    """
    # Track the stages of the synthesis
    epoxide_stage = None
    aminoethanol_coupling_stage = None
    morpholine_formation_stage = None
    morpholine_opening_stage = None

    # Track the presence of morpholine in the final product
    final_product_has_morpholine = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_stage, aminoethanol_coupling_stage, morpholine_formation_stage, morpholine_opening_stage, final_product_has_morpholine

        # Check if this is the final product (depth 0)
        if depth == 0 and node["type"] == "mol":
            final_product_has_morpholine = checker.check_ring("morpholine", node["smiles"])
            print(f"Final product has morpholine: {final_product_has_morpholine}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for epoxide formation
                if checker.check_fg("Oxirane", product) and not any(
                    checker.check_fg("Oxirane", r) for r in reactants
                ):
                    if epoxide_stage is None or depth < epoxide_stage:
                        epoxide_stage = depth
                        print(f"Found epoxide formation at depth {depth}, rsmi: {rsmi}")

                # Check for epoxide presence (as a fallback)
                elif any(checker.check_fg("Oxirane", r) for r in reactants) or checker.check_fg(
                    "Oxirane", product
                ):
                    if epoxide_stage is None or depth < epoxide_stage:
                        epoxide_stage = depth
                        print(f"Found epoxide at depth {depth}, rsmi: {rsmi}")

                # Check for aminoethanol coupling
                # Look for a reaction where both primary amine and primary alcohol are present in reactants
                # and they are being coupled
                if (
                    any(checker.check_fg("Primary alcohol", r) for r in reactants)
                    and any(checker.check_fg("Primary amine", r) for r in reactants)
                    and (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    )
                ):
                    if aminoethanol_coupling_stage is None or depth < aminoethanol_coupling_stage:
                        aminoethanol_coupling_stage = depth
                        print(f"Found aminoethanol coupling at depth {depth}, rsmi: {rsmi}")

                # Check for morpholine formation (ring closure)
                if checker.check_ring("morpholine", product) and not any(
                    checker.check_ring("morpholine", r) for r in reactants
                ):
                    if morpholine_formation_stage is None or depth < morpholine_formation_stage:
                        morpholine_formation_stage = depth
                        print(f"Found morpholine formation at depth {depth}, rsmi: {rsmi}")

                # Check for morpholine ring opening with epoxide
                if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                    if morpholine_opening_stage is None or depth < morpholine_opening_stage:
                        morpholine_opening_stage = depth
                        print(
                            f"Found epoxide ring opening with amine at depth {depth}, rsmi: {rsmi}"
                        )

                # Alternative check for morpholine ring opening
                elif any(
                    checker.check_ring("morpholine", r) for r in reactants
                ) and not checker.check_ring("morpholine", product):
                    if morpholine_opening_stage is None or depth < morpholine_opening_stage:
                        morpholine_opening_stage = depth
                        print(f"Found morpholine ring opening at depth {depth}, rsmi: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if key stages were found
    found_morpholine = morpholine_formation_stage is not None
    found_epoxide = epoxide_stage is not None
    found_opening = morpholine_opening_stage is not None
    found_aminoethanol = aminoethanol_coupling_stage is not None

    print(f"Epoxide stage: {epoxide_stage}")
    print(f"Aminoethanol coupling stage: {aminoethanol_coupling_stage}")
    print(f"Morpholine formation stage: {morpholine_formation_stage}")
    print(f"Morpholine opening stage: {morpholine_opening_stage}")

    # Check for the presence of morpholine in the final product
    if final_product_has_morpholine:
        # If we have a morpholine in the final product and evidence of either:
        # 1. Morpholine formation in the synthesis
        # 2. Epoxide involvement and aminoethanol coupling
        if found_morpholine or (found_epoxide and found_aminoethanol):
            return True
    else:
        # If the final product doesn't have morpholine but we have evidence of
        # morpholine formation and subsequent opening
        if found_morpholine and found_opening:
            return True

    return False
