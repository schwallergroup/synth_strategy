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
    This function detects a synthetic strategy involving sequential formation of multiple pyrazole rings,
    with one formed early in the synthesis and another formed in the final step.
    """
    pyrazole_formations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product from reaction
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a pyrazole formation reaction
            is_pyrazole_reaction = checker.check_reaction("pyrazole", rsmi)

            # Verify by checking pyrazole rings in reactants and product
            reactants_have_pyrazole = any(
                checker.check_ring("pyrazole", smi) for smi in reactants_smiles if smi
            )

            product_has_pyrazole = (
                checker.check_ring("pyrazole", product_smiles) if product_smiles else False
            )

            # If it's explicitly a pyrazole formation reaction or product has pyrazole but reactants don't
            if is_pyrazole_reaction or (product_has_pyrazole and not reactants_have_pyrazole):
                print(f"Pyrazole formation detected at depth {depth}")
                pyrazole_formations.append(depth)

            # Check for hydrazine derivatives in reactants
            has_hydrazine = any(
                checker.check_fg("Hydrazine", smi) for smi in reactants_smiles if smi
            )

            if has_hydrazine:
                print(f"Hydrazine derivative detected in reactants at depth {depth}")

                # Check for carbonyl groups in reactants when hydrazine is present
                has_carbonyl = any(
                    checker.check_fg("Aldehyde", smi) or checker.check_fg("Ketone", smi)
                    for smi in reactants_smiles
                    if smi
                )

                # If we have hydrazine + carbonyl and product has pyrazole, it's likely a pyrazole formation
                if has_carbonyl and product_has_pyrazole and depth not in pyrazole_formations:
                    print(
                        f"Pyrazole formation from hydrazine and carbonyl detected at depth {depth}"
                    )
                    pyrazole_formations.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Pyrazole formations detected at depths: {pyrazole_formations}")

    # Check if we have at least two pyrazole formations at different depths
    if len(pyrazole_formations) >= 2:
        # Lower depths are later stages in synthesis (closer to final product)
        has_final_step_formation = any(depth <= 2 for depth in pyrazole_formations)
        has_early_formation = any(depth > 2 for depth in pyrazole_formations)

        if has_final_step_formation and has_early_formation:
            print("Sequential pyrazole formation strategy detected")
            return True

    return False
