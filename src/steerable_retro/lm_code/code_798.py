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
    This function detects a synthetic strategy involving transformation of an indole ring system.
    This includes both formation and breaking of indole rings.
    """
    indole_transformation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal indole_transformation_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for indole in reactants and product
                indole_in_reactants = any(
                    checker.check_ring("indole", r_smi) for r_smi in reactants_smiles
                )
                indole_in_product = checker.check_ring("indole", product_smiles)

                if indole_in_reactants:
                    print(f"Found indole in reactant(s)")

                if indole_in_product:
                    print(f"Found indole in product")

                # In retrosynthesis, we're looking for:
                # 1. Indole in product but not in reactants (indole formation)
                # 2. Indole in reactants but not in product (indole transformation)

                # Case 1: Indole formation (appears in product but not in reactants)
                if indole_in_product and not indole_in_reactants:
                    print(
                        f"Indole formation detected (indole appears in product but not in reactants)"
                    )

                    # Check if this is a known indole-forming reaction
                    indole_forming_reactions = [
                        "Fischer indole",
                        "{Fischer indole}",
                        "Sonogashira",
                        "{indole}",
                    ]
                    if any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in indole_forming_reactions
                    ):
                        print(f"Detected known indole-forming reaction")
                        indole_transformation_detected = True

                    # Check for cyclization reactions that might form indole
                    cyclization_reactions = [
                        "Sonogashira acetylene_aryl halide",
                        "Sonogashira alkyne_aryl halide",
                    ]
                    if any(
                        checker.check_reaction(rxn_type, rsmi) for rxn_type in cyclization_reactions
                    ):
                        print(f"Detected cyclization reaction that could form indole")
                        indole_transformation_detected = True

                    # Check for precursors that might form indole
                    alkyne_in_reactants = any(
                        checker.check_fg("Alkyne", r_smi) for r_smi in reactants_smiles
                    )
                    aniline_in_reactants = any(
                        checker.check_fg("Aniline", r_smi) for r_smi in reactants_smiles
                    )

                    if alkyne_in_reactants and aniline_in_reactants:
                        print(f"Alkyne and aniline found in reactants, potential indole formation")
                        indole_transformation_detected = True
                    elif aniline_in_reactants:
                        print(f"Aniline found in reactants, potential indole formation")
                        indole_transformation_detected = True
                    elif alkyne_in_reactants:
                        print(f"Alkyne found in reactants, potential indole formation with amine")
                        # Check if any reactant has an amine group
                        amine_in_reactants = any(
                            any(
                                checker.check_fg(fg, r_smi)
                                for fg in ["Primary amine", "Secondary amine"]
                            )
                            for r_smi in reactants_smiles
                        )
                        if amine_in_reactants:
                            print(f"Amine and alkyne found, potential indole formation")
                            indole_transformation_detected = True

                # Case 2: Indole transformation (disappears from reactants to product)
                elif indole_in_reactants and not indole_in_product:
                    print(
                        f"Indole transformation detected (indole disappears from reactants to product)"
                    )

                    # Check if product has aniline or amine groups (common after indole breaking)
                    if checker.check_fg("Aniline", product_smiles):
                        print(f"Aniline found in product after indole transformation")
                        indole_transformation_detected = True
                    elif any(
                        checker.check_fg(fg, product_smiles)
                        for fg in ["Primary amine", "Secondary amine", "Tertiary amine"]
                    ):
                        print(f"Amine group found in product after indole transformation")
                        indole_transformation_detected = True

                    # Check if this is a known ring-breaking reaction
                    ring_breaking_reactions = [
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Reductive amination",
                    ]
                    if any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in ring_breaking_reactions
                    ):
                        print(f"Detected known ring-breaking reaction")
                        indole_transformation_detected = True

                # Case 3: Indole modification (present in both but modified)
                elif indole_in_reactants and indole_in_product:
                    # Check for functional group changes on the indole
                    indole_modifying_reactions = [
                        "Aromatic nitration with HNO3",
                        "Aromatic nitration with NO3 salt",
                        "Aromatic nitration with NO2 salt",
                        "Aromatic nitration with alkyl NO2",
                        "Reduction of nitro groups to amines",
                        "Esterification of Carboxylic Acids",
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        "Oxidation of aldehydes to carboxylic acids",
                        "Reduction of nitrile to amine",
                    ]

                    if any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in indole_modifying_reactions
                    ):
                        print(f"Detected reaction that modifies indole structure")
                        indole_transformation_detected = True

                    # Check for specific functional group changes on indole
                    functional_groups = [
                        "Nitro group",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Carboxylic acid",
                        "Ester",
                        "Aldehyde",
                        "Ketone",
                        "Nitrile",
                        "Aromatic halide",
                        "Primary alcohol",
                        "Secondary alcohol",
                    ]

                    for fg in functional_groups:
                        fg_in_product = checker.check_fg(fg, product_smiles)
                        fg_in_reactants = any(
                            checker.check_fg(fg, r_smi) for r_smi in reactants_smiles
                        )
                        if fg_in_product != fg_in_reactants:
                            print(f"Functional group change detected: {fg}")
                            indole_transformation_detected = True
                            break

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Indole transformation detected: {indole_transformation_detected}")
    return indole_transformation_detected
