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
    This function detects a synthetic strategy involving late-stage N-oxide reduction
    in a pyridine-containing substrate.
    """
    # Initialize tracking variables
    has_noxide_reduction = False
    depth_of_noxide_reduction = float("inf")

    def is_pyridine_noxide(smiles):
        """Check if the molecule contains a pyridine N-oxide structure"""
        # Check for pyridine first
        if not checker.check_ring("pyridine", smiles):
            return False

        # Check for N-oxide pattern in SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Look for [n+]([O-]) pattern which represents N-oxide in pyridine
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "N" and atom.GetFormalCharge() == 1:
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == "O" and neighbor.GetFormalCharge() == -1:
                            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_noxide_reduction, depth_of_noxide_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # In retrosynthesis, the product of the forward reaction is the starting material
                # and the reactants are what we're synthesizing from

                # Check if product has pyridine N-oxide (forward direction)
                product_has_pyridine_noxide = is_pyridine_noxide(product_smiles)

                # Check if any reactant has pyridine but not N-oxide
                reactant_has_pyridine_without_noxide = False
                for reactant_smiles in reactants_smiles:
                    if checker.check_ring("pyridine", reactant_smiles) and not is_pyridine_noxide(
                        reactant_smiles
                    ):
                        reactant_has_pyridine_without_noxide = True
                        break

                print(f"Product has pyridine N-oxide: {product_has_pyridine_noxide}")
                print(
                    f"Reactant has pyridine without N-oxide: {reactant_has_pyridine_without_noxide}"
                )

                # In forward direction: pyridine → pyridine N-oxide is oxidation
                # In retrosynthesis: we're looking for pyridine N-oxide → pyridine (reduction)
                if product_has_pyridine_noxide and reactant_has_pyridine_without_noxide:
                    print("Detected N-oxide formation (oxidation) in forward direction")
                    print("This corresponds to N-oxide reduction in retrosynthesis")

                    # This is a late-stage N-oxide reduction in retrosynthetic analysis
                    has_noxide_reduction = True
                    depth_of_noxide_reduction = min(depth_of_noxide_reduction, depth)
                    print(f"Confirmed N-oxide reduction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a late-stage reduction (depth 0 or 1)
    is_late_stage = has_noxide_reduction and depth_of_noxide_reduction <= 1

    if is_late_stage:
        print(f"Detected late-stage N-oxide reduction at depth {depth_of_noxide_reduction}")
    else:
        if has_noxide_reduction:
            print(f"N-oxide reduction found but not late-stage (depth {depth_of_noxide_reduction})")
        else:
            print("No N-oxide reduction found in the route")

    return is_late_stage
