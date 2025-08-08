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
    This function detects a synthetic strategy involving a sequence of
    transformations from nitrile to hydroxylamine to heterocycle.
    """
    # Track the sequence of transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrile to hydroxylamine/amidoxime transformation
                if checker.check_fg("Nitrile", reactants_smiles) and not checker.check_fg(
                    "Nitrile", product_smiles
                ):
                    # Check for specific reactions
                    if checker.check_reaction(
                        "N-hydroxyimidamide from nitrile and hydroxylamine", rsmi
                    ) or checker.check_reaction("Amidoxime from nitrile and hydroxylamine", rsmi):
                        print(
                            f"Detected specific nitrile to hydroxylamine reaction at depth {depth}"
                        )
                        print(f"Reactants: {reactants_smiles}")
                        print(f"Product: {product_smiles}")
                        transformations.append(("nitrile_to_hydroxylamine", depth))
                    # Check for hydroxylamine or amidoxime formation using patterns
                    elif (
                        checker.check_fg("Oxime", product_smiles)
                        or "NOH" in product_smiles
                        or "N(O)" in product_smiles
                        or "N-O" in product_smiles
                    ):
                        print(f"Potential nitrile to hydroxylamine transformation at depth {depth}")
                        print(f"Reactants: {reactants_smiles}")
                        print(f"Product: {product_smiles}")
                        transformations.append(("nitrile_to_hydroxylamine", depth))

                # Check for hydroxylamine to heterocycle transformation
                if (
                    checker.check_fg("Oxime", reactants_smiles)
                    or "NOH" in reactants_smiles
                    or "N(O)" in reactants_smiles
                    or "N-O" in reactants_smiles
                ):

                    # Check for specific heterocycle formation reactions
                    if (
                        checker.check_reaction(
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                        or checker.check_reaction(
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                        )
                        or checker.check_reaction(
                            "Azide-nitrile click cycloaddition to tetrazole", rsmi
                        )
                        or checker.check_reaction(
                            "Azide-nitrile click cycloaddition to triazole", rsmi
                        )
                        or checker.check_reaction("Pyrazole formation", rsmi)
                        or checker.check_reaction(
                            "1,2,4-oxadiazol-5(2H)-one synthesis from nitrile, hydrogen carbonate, and hydroxylamine",
                            rsmi,
                        )
                    ):

                        print(
                            f"Detected specific hydroxylamine to heterocycle reaction at depth {depth}"
                        )
                        print(f"Reactants: {reactants_smiles}")
                        print(f"Product: {product_smiles}")
                        transformations.append(
                            ("hydroxylamine_to_heterocycle", depth, "specific_reaction")
                        )

                    # Check for various heterocycles in the product
                    heterocycle_found = False
                    for ring_name in [
                        "isoxazole",
                        "oxadiazole",
                        "tetrazole",
                        "triazole",
                        "oxazole",
                        "pyrazole",
                    ]:
                        if checker.check_ring(ring_name, product_smiles) and not checker.check_ring(
                            ring_name, reactants_smiles
                        ):
                            print(f"Hydroxylamine to {ring_name} transformation at depth {depth}")
                            print(f"Reactants: {reactants_smiles}")
                            print(f"Product: {product_smiles}")
                            transformations.append(
                                ("hydroxylamine_to_heterocycle", depth, ring_name)
                            )
                            heterocycle_found = True
                            break

                    # If no specific ring was found but a cyclic structure appears
                    if not heterocycle_found and any(c in product_smiles for c in "123456789"):
                        print(
                            f"Potential hydroxylamine to heterocycle transformation at depth {depth}"
                        )
                        print(f"Reactants: {reactants_smiles}")
                        print(f"Product: {product_smiles}")
                        transformations.append(("hydroxylamine_to_heterocycle", depth, "unknown"))

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have the complete sequence
    nitrile_to_hydroxylamine_steps = [
        (t[0], t[1]) for t in transformations if t[0] == "nitrile_to_hydroxylamine"
    ]
    hydroxylamine_to_heterocycle_steps = [
        (t[0], t[1]) for t in transformations if t[0] == "hydroxylamine_to_heterocycle"
    ]

    # Check if both transformations are present
    if nitrile_to_hydroxylamine_steps and hydroxylamine_to_heterocycle_steps:
        # Check if the sequence is in the correct order (nitrile → hydroxylamine → heterocycle)
        # In retrosynthetic analysis, higher depth means earlier in the synthesis
        nitrile_depths = [t[1] for t in nitrile_to_hydroxylamine_steps]
        heterocycle_depths = [t[1] for t in hydroxylamine_to_heterocycle_steps]

        # Check if at least one nitrile step occurs before a heterocycle step
        sequence_in_order = any(
            nitrile_depth > heterocycle_depth
            for nitrile_depth in nitrile_depths
            for heterocycle_depth in heterocycle_depths
        )

        strategy_present = sequence_in_order
        print(f"Nitrile to heterocycle sequence detected: {strategy_present}")
        print(f"Nitrile steps at depths: {nitrile_depths}")
        print(f"Heterocycle steps at depths: {heterocycle_depths}")
        return strategy_present
    else:
        print(f"Nitrile to heterocycle sequence detected: False")
        print(
            f"Missing steps: {'nitrile_to_hydroxylamine' if not nitrile_to_hydroxylamine_steps else ''} {'hydroxylamine_to_heterocycle' if not hydroxylamine_to_heterocycle_steps else ''}"
        )
        return False
