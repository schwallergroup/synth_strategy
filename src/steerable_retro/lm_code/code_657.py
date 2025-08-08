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
    This function detects if the synthesis involves building a quinazoline scaffold
    followed by functional group modifications and late-stage coupling.
    """
    # Track the key features of this strategy
    has_quinazoline_formation = False
    has_functional_group_modification = False
    has_late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_quinazoline_formation, has_functional_group_modification, has_late_stage_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for quinazoline formation (early stage)
                if depth >= 2:  # Early stage
                    # Check if product has quinazoline
                    if checker.check_ring("quinoline", product_smiles):
                        # Check if reactants don't have quinazoline
                        reactants_have_quinazoline = False
                        for reactant_smiles in reactants_smiles:
                            if checker.check_ring("quinoline", reactant_smiles):
                                reactants_have_quinazoline = True
                                break

                        if not reactants_have_quinazoline:
                            print(f"Quinazoline formation detected at depth {depth}")
                            # Check if it's specifically a Niementowski quinazoline formation
                            if checker.check_reaction("{Niementowski_quinazoline}", rsmi):
                                print(
                                    f"Niementowski quinazoline formation confirmed at depth {depth}"
                                )
                                has_quinazoline_formation = True
                            else:
                                # Alternative check for quinazoline formation
                                has_quinazoline_formation = True
                                print(f"General quinazoline formation detected at depth {depth}")

                # Check for functional group modifications (middle stage)
                if 1 <= depth < 3:  # Middle stage
                    # Check for alcohol to halide conversion
                    if (
                        any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            for r in reactants_smiles
                        )
                        and checker.check_fg("Primary halide", product_smiles)
                        or checker.check_fg("Secondary halide", product_smiles)
                        or checker.check_fg("Tertiary halide", product_smiles)
                    ):

                        # Verify it's a specific reaction type
                        if (
                            checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                            or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                            or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                            or checker.check_reaction("Appel reaction", rsmi)
                        ):
                            print(f"Alcohol to halide conversion detected at depth {depth}")
                            has_functional_group_modification = True

                    # Check for nitrile to amide conversion
                    if checker.check_reaction("Nitrile to amide", rsmi):
                        print(f"Nitrile to amide conversion detected at depth {depth}")
                        has_functional_group_modification = True

                    # Check for other common functional group modifications
                    if (
                        checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                        or checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                        or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    ):
                        print(f"Other functional group modification detected at depth {depth}")
                        has_functional_group_modification = True

                # Check for late-stage coupling
                if depth <= 1:  # Late stage
                    # Check for specific coupling reactions
                    if (
                        checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction("{Buchwald-Hartwig}", rsmi)
                        or checker.check_reaction("Goldberg coupling", rsmi)
                    ):
                        print(f"Late-stage aromatic amine coupling detected at depth {depth}")
                        has_late_stage_coupling = True

                    # Alternative check for coupling reactions
                    if len(reactants_smiles) >= 2:
                        # Check if both reactants have aromatic rings
                        aromatic_reactants = 0
                        for reactant_smiles in reactants_smiles:
                            if any(
                                checker.check_ring(ring, reactant_smiles)
                                for ring in ["benzene", "pyridine", "quinoline"]
                            ):
                                aromatic_reactants += 1

                        # If we have at least 2 aromatic reactants and the product has an aromatic amine
                        if aromatic_reactants >= 2:
                            if checker.check_fg("Aniline", product_smiles) or checker.check_fg(
                                "Secondary amine", product_smiles
                            ):
                                print(f"Late-stage aromatic coupling detected at depth {depth}")
                                has_late_stage_coupling = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is present if all three key features are detected
    print(f"Quinazoline formation: {has_quinazoline_formation}")
    print(f"Functional group modification: {has_functional_group_modification}")
    print(f"Late-stage coupling: {has_late_stage_coupling}")

    # Check if any two of the three features are detected (relaxed condition)
    features_count = sum(
        [has_quinazoline_formation, has_functional_group_modification, has_late_stage_coupling]
    )
    return features_count >= 2  # Return true if at least 2 features are detected
