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
    Detects a synthetic strategy involving thiazole formation from α-bromoketone and thiourea,
    with the α-bromoketone being formed from a ketone, which in turn comes from a Weinreb amide.
    """
    # Initialize tracking variables
    has_thiazole_formation = False
    has_alpha_bromination = False
    has_weinreb_amide = False

    def dfs_traverse(node):
        nonlocal has_thiazole_formation, has_alpha_bromination, has_weinreb_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for thiazole formation
                product_has_thiazole = checker.check_ring("thiazole", product_smiles)
                reactants_have_thiazole = any(
                    checker.check_ring("thiazole", r) for r in reactants_smiles
                )
                reactants_have_thiourea = any(
                    checker.check_fg("Thiourea", r) for r in reactants_smiles
                )

                # Check for halide or alpha-haloketone in reactants
                reactants_have_halide = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Aromatic halide", r)
                    for r in reactants_smiles
                )

                reactants_have_ketone = any(checker.check_fg("Ketone", r) for r in reactants_smiles)

                # Check for alpha-haloketone specifically
                alpha_haloketone_present = False
                for r in reactants_smiles:
                    if checker.check_fg("Ketone", r) and (
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Aromatic halide", r)
                    ):
                        alpha_haloketone_present = True
                        break

                if (
                    product_has_thiazole
                    and not reactants_have_thiazole
                    and reactants_have_thiourea
                    and (
                        alpha_haloketone_present
                        or (reactants_have_halide and reactants_have_ketone)
                    )
                ):
                    print("Detected thiazole formation from thiourea and α-bromoketone")
                    has_thiazole_formation = True

                # Check for α-bromination of ketone
                if (
                    checker.check_reaction("Wohl-Ziegler bromination carbonyl primary", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination carbonyl secondary", rsmi)
                    or checker.check_reaction("Bromination", rsmi)
                    or checker.check_reaction("Aromatic bromination", rsmi)
                ):

                    # Verify reactants contain ketone and products contain halogenated compound
                    reactant_has_ketone = any(
                        checker.check_fg("Ketone", r) for r in reactants_smiles
                    )
                    product_has_halide = (
                        checker.check_fg("Primary halide", product_smiles)
                        or checker.check_fg("Secondary halide", product_smiles)
                        or checker.check_fg("Aromatic halide", product_smiles)
                    )

                    if reactant_has_ketone and product_has_halide:
                        print("Detected α-bromination of ketone")
                        has_alpha_bromination = True

                # Additional check for bromination without specific reaction type
                if not has_alpha_bromination:
                    reactant_has_ketone = any(
                        checker.check_fg("Ketone", r) for r in reactants_smiles
                    )
                    product_has_ketone = checker.check_fg("Ketone", product_smiles)
                    product_has_halide = checker.check_fg(
                        "Primary halide", product_smiles
                    ) or checker.check_fg("Secondary halide", product_smiles)

                    if reactant_has_ketone and product_has_ketone and product_has_halide:
                        print("Detected α-bromination of ketone (alternative check)")
                        has_alpha_bromination = True

                # Check for Weinreb amide conversion
                if checker.check_reaction("Ketone from Weinreb amide", rsmi):
                    print("Detected Weinreb amide conversion to ketone")
                    has_weinreb_amide = True

                # Alternative check for Weinreb amide conversion if reaction type check fails
                if not has_weinreb_amide:
                    # Check for Weinreb amide structure in reactants
                    reactants_have_weinreb = any(
                        "N(C)OC" in r or "N(OC)" in r or checker.check_fg("Tertiary amide", r)
                        for r in reactants_smiles
                    )
                    product_has_ketone = checker.check_fg("Ketone", product_smiles)

                    if reactants_have_weinreb and product_has_ketone:
                        print("Detected Weinreb amide conversion to ketone (alternative check)")
                        has_weinreb_amide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the complete strategy is present
    strategy_detected = has_thiazole_formation and has_alpha_bromination and has_weinreb_amide
    if strategy_detected:
        print("Complete thiazole from bromoketone-thiourea strategy detected")
    else:
        print(
            f"Strategy incomplete: thiazole formation: {has_thiazole_formation}, α-bromination: {has_alpha_bromination}, Weinreb amide: {has_weinreb_amide}"
        )

    return strategy_detected
