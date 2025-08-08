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
    Detects a strategy involving pyrazole formation from a diazonium intermediate.
    """
    # Track where diazonium/diazo groups and pyrazole formations occur
    diazonium_depths = []
    pyrazole_formation_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for diazo/diazonium group in reactants
                diazo_present = False
                for reactant in reactants:
                    if (
                        checker.check_fg("Diazo", reactant)
                        or checker.check_fg("Diazene", reactant)
                        or "N#N" in reactant
                    ):
                        diazo_present = True
                        diazonium_depths.append(depth)
                        print(f"Detected diazo/diazonium group at depth {depth} in: {reactant}")
                        break

                # Check for alkyne in reactants
                alkyne_present = any(checker.check_fg("Alkyne", reactant) for reactant in reactants)
                if alkyne_present:
                    print(f"Detected alkyne at depth {depth}")

                # Check for pyrazole in product
                pyrazole_in_product = checker.check_ring("pyrazole", product)
                if pyrazole_in_product:
                    print(f"Detected pyrazole ring in product at depth {depth}")

                # Check for various cycloaddition reactions that could form pyrazoles
                is_cycloaddition = (
                    checker.check_reaction("Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of diazoalkane and alkyne", rsmi)
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkyne", rsmi)
                    or checker.check_reaction("pyrazole", rsmi)
                    or checker.check_reaction("{pyrazole}", rsmi)
                )

                if is_cycloaddition:
                    print(f"Detected cycloaddition reaction at depth {depth}")

                # Check if this reaction forms a pyrazole
                # Either through a specific cycloaddition or if we have both diazo and pyrazole
                if (alkyne_present and pyrazole_in_product and is_cycloaddition) or (
                    diazo_present and pyrazole_in_product
                ):
                    pyrazole_formation_depths.append(depth)
                    print(f"Detected pyrazole formation at depth {depth}")

                # Also check if we have a diazo compound reacting with an alkyne to form any product
                # This catches cases where the reaction might not be explicitly labeled
                if diazo_present and alkyne_present and not is_cycloaddition:
                    # Check if product contains a newly formed pyrazole
                    if pyrazole_in_product and not any(
                        checker.check_ring("pyrazole", reactant) for reactant in reactants
                    ):
                        pyrazole_formation_depths.append(depth)
                        print(
                            f"Detected implicit pyrazole formation from diazo and alkyne at depth {depth}"
                        )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both diazonium/diazo groups and pyrazole formation in the route
    if diazonium_depths and pyrazole_formation_depths:
        # In retrosynthetic analysis, we need to check if the diazonium appears at the same depth
        # or deeper (earlier in synthesis) than the pyrazole formation
        strategy_detected = False
        for d_depth in diazonium_depths:
            for p_depth in pyrazole_formation_depths:
                # If diazonium is at same depth or deeper than pyrazole formation,
                # it's part of the strategy (in retrosynthesis, deeper = earlier in forward synthesis)
                if d_depth >= p_depth:
                    strategy_detected = True
                    print(
                        f"Confirmed diazonium-pyrazole formation strategy: diazonium at depth {d_depth}, pyrazole formation at depth {p_depth}"
                    )

        if strategy_detected:
            print("Detected diazonium-pyrazole formation strategy")
            return True

    if diazonium_depths:
        print("Found diazonium/diazo group but no complete pyrazole formation strategy")
    else:
        print("No diazonium/diazo groups found in the route")

    return False
