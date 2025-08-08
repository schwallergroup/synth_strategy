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
    This function detects a synthetic strategy involving late-stage
    conversion of a nitrile to a ketone.
    """
    # Initialize flag
    has_late_stage_nitrile_to_ketone = False

    # Track nitrile-containing molecules and their depths
    nitrile_molecules = {}  # {smiles: depth}

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_nitrile_to_ketone

        # For molecule nodes, check and track nitrile-containing molecules
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Nitrile", mol_smiles):
                nitrile_molecules[mol_smiles] = depth
                print(f"Found nitrile-containing molecule at depth {depth}: {mol_smiles}")

        # For reaction nodes, check for nitrile to ketone transformation
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth ≤ 3)
            if depth <= 3:
                # Check if any reactant contains nitrile group
                reactants = reactants_str.split(".")
                nitrile_reactants = [r for r in reactants if checker.check_fg("Nitrile", r)]

                # Check if product contains ketone group
                has_ketone_product = checker.check_fg("Ketone", product_str)

                # Direct nitrile to ketone reactions
                if nitrile_reactants and has_ketone_product:
                    # Check for Grignard reaction (most common nitrile to ketone)
                    if checker.check_reaction("Grignard from nitrile to ketone", rsmi):
                        print(
                            f"Found late-stage direct nitrile to ketone transformation at depth {depth}: {rsmi}"
                        )
                        print(f"  Reaction type: Grignard from nitrile to ketone")
                        print(f"  Nitrile reactants: {nitrile_reactants}")
                        print(f"  Ketone product: {product_str}")
                        has_late_stage_nitrile_to_ketone = True

                    # Check for hydrolysis/hydration of nitrile to ketone
                    elif any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Nef reaction (nitro to ketone)",  # Some nitriles can be converted via nitro compounds
                            "Hydration of alkyne to ketone",  # Some nitriles can be converted via alkynes
                        ]
                    ):
                        rxn_type = next(
                            rxn
                            for rxn in [
                                "Nef reaction (nitro to ketone)",
                                "Hydration of alkyne to ketone",
                            ]
                            if checker.check_reaction(rxn, rsmi)
                        )
                        print(
                            f"Found late-stage nitrile to ketone transformation at depth {depth}: {rsmi}"
                        )
                        print(f"  Reaction type: {rxn_type}")
                        print(f"  Nitrile reactants: {nitrile_reactants}")
                        print(f"  Ketone product: {product_str}")
                        has_late_stage_nitrile_to_ketone = True

                    # Check for other potential nitrile to ketone conversions
                    elif not any(
                        checker.check_fg(fg, product_str)
                        for fg in ["Carboxylic acid", "Ester", "Amide"]
                    ):
                        # If product doesn't contain other common nitrile derivatives, likely a ketone conversion
                        print(
                            f"Found potential late-stage nitrile to ketone transformation at depth {depth}: {rsmi}"
                        )
                        print(f"  Nitrile reactants: {nitrile_reactants}")
                        print(f"  Ketone product: {product_str}")
                        has_late_stage_nitrile_to_ketone = True

                # Check for two-step nitrile to ketone via amide or carboxylic acid
                elif not nitrile_reactants and has_ketone_product:
                    # Check if reactants contain amide or carboxylic acid (potential nitrile derivatives)
                    amide_reactants = [
                        r
                        for r in reactants
                        if checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                    ]

                    acid_reactants = [
                        r for r in reactants if checker.check_fg("Carboxylic acid", r)
                    ]

                    # Check for amide to ketone conversion
                    if amide_reactants and any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Ketone from Weinreb amide",
                            "Asymmetric ketones from N,N-dimethylamides",
                        ]
                    ):
                        rxn_type = next(
                            rxn
                            for rxn in [
                                "Ketone from Weinreb amide",
                                "Asymmetric ketones from N,N-dimethylamides",
                            ]
                            if checker.check_reaction(rxn, rsmi)
                        )

                        # Check if the amide was derived from a nitrile in a previous step
                        for amide in amide_reactants:
                            for nitrile_mol, nitrile_depth in nitrile_molecules.items():
                                if nitrile_depth < depth and depth <= 3:
                                    print(
                                        f"Found late-stage indirect nitrile to ketone via amide at depth {depth}: {rsmi}"
                                    )
                                    print(f"  Reaction type: {rxn_type}")
                                    print(f"  Amide reactant (from nitrile): {amide}")
                                    print(f"  Ketone product: {product_str}")
                                    has_late_stage_nitrile_to_ketone = True

                    # Check for acid to ketone conversion
                    elif acid_reactants and any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Ketonization by decarboxylation of carbonic acids",
                            "Ketonization by decarboxylation of acid halides",
                        ]
                    ):
                        rxn_type = next(
                            rxn
                            for rxn in [
                                "Ketonization by decarboxylation of carbonic acids",
                                "Ketonization by decarboxylation of acid halides",
                            ]
                            if checker.check_reaction(rxn, rsmi)
                        )

                        # Check if the acid was derived from a nitrile in a previous step
                        for acid in acid_reactants:
                            for nitrile_mol, nitrile_depth in nitrile_molecules.items():
                                if nitrile_depth < depth and depth <= 3:
                                    print(
                                        f"Found late-stage indirect nitrile to ketone via acid at depth {depth}: {rsmi}"
                                    )
                                    print(f"  Reaction type: {rxn_type}")
                                    print(f"  Acid reactant (from nitrile): {acid}")
                                    print(f"  Ketone product: {product_str}")
                                    has_late_stage_nitrile_to_ketone = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Strategy detection results:")
    print(f"  Late-stage nitrile to ketone: {has_late_stage_nitrile_to_ketone}")

    return has_late_stage_nitrile_to_ketone
