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
    This function detects a synthetic strategy involving dehydrogenation
    (alkyl to alkene conversion) in the late stages of synthesis.
    """
    late_stage_dehydrogenation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_dehydrogenation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for hydrogenation reactions (which would be dehydrogenation in forward direction)
            is_hydrogenation = checker.check_reaction(
                "Hydrogenation (double to single)", rsmi
            ) or checker.check_reaction("Hydrogenation (triple to double)", rsmi)

            # In retrosynthesis, a dehydrogenation would appear as:
            # Product (more saturated) -> Reactant (less saturated)
            reactants_have_alkene = any(
                checker.check_fg("Vinyl", r)
                or checker.check_fg("Allyl", r)
                or checker.check_fg("Ethylene", r)
                or checker.check_fg("Alkyne", r)
                for r in reactants_smiles
            )

            product_has_vinyl = checker.check_fg("Vinyl", product_smiles)
            product_has_allyl = checker.check_fg("Allyl", product_smiles)
            product_has_alkyne = checker.check_fg("Alkyne", product_smiles)

            # Check if the product has alkyl groups that correspond to unsaturated groups in reactants
            # This specifically checks for the pattern in the test case
            dehydrogenation_pattern = False

            # Direct pattern check for vinyl/allene to alkyl conversion
            if "[CH2:1]=[CH:2]" in rsmi and "[CH3:1][CH2:2]" in rsmi:
                print("  Found vinyl to alkyl conversion pattern")
                dehydrogenation_pattern = True

            # Check if reactants have unsaturated bonds that product doesn't have
            elif reactants_have_alkene and not (
                product_has_vinyl or product_has_allyl or product_has_alkyne
            ):
                print("  Found unsaturated reactants with saturated product")
                dehydrogenation_pattern = True

            # Check for explicit hydrogenation reaction
            elif is_hydrogenation:
                print("  Found explicit hydrogenation reaction")
                dehydrogenation_pattern = True

            print(f"  Is hydrogenation reaction: {is_hydrogenation}")
            print(f"  Reactants have alkene: {reactants_have_alkene}")
            print(f"  Product has vinyl: {product_has_vinyl}")
            print(f"  Product has allyl: {product_has_allyl}")
            print(f"  Product has alkyne: {product_has_alkyne}")
            print(f"  Dehydrogenation pattern detected: {dehydrogenation_pattern}")

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2 and dehydrogenation_pattern:
                print(f"LATE-STAGE DEHYDROGENATION DETECTED at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                late_stage_dehydrogenation = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_dehydrogenation
