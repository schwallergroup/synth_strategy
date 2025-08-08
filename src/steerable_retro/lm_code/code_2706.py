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
    Detects a strategy where a heterocyclic scaffold is constructed early in the synthesis,
    followed by sequential nucleophilic aromatic substitutions of halogens with nitrogen nucleophiles.
    """
    # Track key features
    heterocycle_formation = False
    heterocycle_depth = -1
    halogen_displacements = 0
    halogen_displacement_depths = []

    # List of heterocyclic rings to check
    heterocycle_rings = [
        "quinoline",
        "isoquinoline",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "purine",
        "pyrrole",
        "pyrazole",
        "furan",
        "thiophene",
        "oxadiazole",
        "thiadiazole",
    ]

    def is_heterocycle_formation(rsmi):
        """Check if the reaction forms a heterocycle"""
        reactants_smiles = rsmi.split(">")[0].split(".")
        product_smiles = rsmi.split(">")[-1]

        # Check if product contains a heterocycle that wasn't in any reactant
        for ring in heterocycle_rings:
            if checker.check_ring(ring, product_smiles):
                # Check if any reactant already had this ring
                if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                    print(f"Detected heterocycle formation: {ring}")
                    return True

        # Check for specific heterocycle formation reactions
        heterocycle_forming_reactions = [
            "Friedlaender chinoline",
            "Fischer indole",
            "benzimidazole_derivatives_carboxylic-acid/ester",
            "benzimidazole_derivatives_aldehyde",
            "benzothiazole",
            "benzoxazole_arom-aldehyde",
            "benzoxazole_carboxylic-acid",
            "thiazole",
            "Paal-Knorr pyrrole",
            "Formation of NOS Heterocycles",
            "Pictet-Spengler",
            "benzimidazole",
            "benzoxazole",
            "indole",
            "oxadiazole",
            "pyrazole",
            "tetrazole_terminal",
            "tetrazole_connect_regioisomere_1",
            "tetrazole_connect_regioisomere_2",
            "1,2,4-triazole_acetohydrazide",
            "1,2,4-triazole_carboxylic-acid/ester",
            "3-nitrile-pyridine",
            "triaryl-imidazole",
            "benzofuran",
            "benzothiophene",
        ]

        for rxn_type in heterocycle_forming_reactions:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected heterocycle formation reaction: {rxn_type}")
                return True

        return False

    def is_halogen_displacement(rsmi):
        """Check if the reaction is a halogen displacement with a nitrogen nucleophile"""
        # Check for specific nucleophilic aromatic substitution reactions
        displacement_reactions = [
            "Goldberg coupling",
            "Ullmann-Goldberg Substitution amine",
            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Buchwald-Hartwig",
            "heteroaromatic_nuc_sub",
            "nucl_sub_aromatic_ortho_nitro",
            "nucl_sub_aromatic_para_nitro",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
            "Goldberg coupling aryl amine-aryl chloride",
            "Goldberg coupling aryl amide-aryl chloride",
            "Ullmann condensation",
            "N-arylation_heterocycles",
        ]

        for rxn_type in displacement_reactions:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected halogen displacement reaction: {rxn_type}")
                return True

        # If no specific reaction type matches, check for the pattern manually
        reactants_smiles = rsmi.split(">")[0].split(".")
        product_smiles = rsmi.split(">")[-1]

        # Check if any reactant has an aromatic halide
        has_aromatic_halide = False
        for r in reactants_smiles:
            if checker.check_fg("Aromatic halide", r):
                has_aromatic_halide = True
                break

        # Check if any reactant has a nitrogen nucleophile
        has_nitrogen_nucleophile = False
        for r in reactants_smiles:
            if (
                checker.check_fg("Primary amine", r)
                or checker.check_fg("Secondary amine", r)
                or checker.check_fg("Tertiary amine", r)
                or checker.check_fg("Aniline", r)
                or checker.check_fg("Azide", r)
            ):
                has_nitrogen_nucleophile = True
                break

        if has_aromatic_halide and has_nitrogen_nucleophile:
            # Check if product has a new C-N bond (aniline or similar)
            if (
                checker.check_fg("Aniline", product_smiles)
                or checker.check_fg("Secondary amine", product_smiles)
                or checker.check_fg("Tertiary amine", product_smiles)
            ):
                # Verify that a halogen was actually displaced
                if not checker.check_fg("Aromatic halide", product_smiles) or sum(
                    1 for r in reactants_smiles if checker.check_fg("Aromatic halide", r)
                ) > (1 if checker.check_fg("Aromatic halide", product_smiles) else 0):
                    print("Detected halogen displacement pattern")
                    return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation, heterocycle_depth, halogen_displacements, halogen_displacement_depths

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for heterocycle formation
                if not heterocycle_formation and is_heterocycle_formation(rsmi):
                    heterocycle_formation = True
                    heterocycle_depth = depth
                    print(f"Detected heterocycle formation at depth {depth}")

                # Check for halogen displacement with nitrogen nucleophile
                if is_halogen_displacement(rsmi):
                    halogen_displacements += 1
                    halogen_displacement_depths.append(depth)
                    print(f"Detected halogen displacement at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need heterocycle formation early (higher depth) and at least 1 halogen displacement later
    strategy_present = False

    print(f"Heterocycle formation: {heterocycle_formation} at depth {heterocycle_depth}")
    print(f"Halogen displacements: {halogen_displacements} at depths {halogen_displacement_depths}")

    if heterocycle_formation and halogen_displacements >= 1:
        # In retrosynthesis, early-stage reactions (heterocycle formation) should have higher depth
        # than late-stage reactions (halogen displacements)
        if all(heterocycle_depth > d for d in halogen_displacement_depths):
            print("Heterocycle formation occurs earlier in synthesis than halogen displacements")

            # Check if there are at least two halogen displacements or one significant one
            if halogen_displacements >= 2:
                # Sort depths in descending order for retrosynthesis (higher depth = earlier in synthesis)
                sorted_depths = sorted(halogen_displacement_depths, reverse=True)

                # Check if there are sequential displacements (allowing for intermediate steps)
                for i in range(len(sorted_depths) - 1):
                    if abs(sorted_depths[i] - sorted_depths[i + 1]) <= 3:  # Allow more flexibility
                        strategy_present = True
                        print(
                            f"Found sequential halogen displacements at depths {sorted_depths[i]} and {sorted_depths[i+1]}"
                        )
                        break
            else:
                # If there's only one displacement but it's significant (late in synthesis)
                # and the heterocycle was formed early, still consider the strategy present
                strategy_present = True
                print(
                    "Only one halogen displacement detected, but considering strategy present due to clear pattern"
                )

    print(f"Strategy present: {strategy_present}")
    return strategy_present
