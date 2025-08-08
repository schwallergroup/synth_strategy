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
    This function detects if the synthetic route maintains a cyclopropyl group throughout the synthesis.
    """
    # Track if cyclopropyl is present in the target molecule
    target_has_cyclopropyl = False

    # Track if the cyclopropyl group is maintained throughout the synthesis
    cyclopropyl_maintained = True

    # Track molecules that should have cyclopropyl based on synthetic path
    molecules_with_cyclopropyl = set()

    # Track reactions where cyclopropyl is created
    cyclopropyl_creation_reactions = set()

    def dfs_traverse(node, depth=0):
        nonlocal target_has_cyclopropyl, cyclopropyl_maintained, molecules_with_cyclopropyl

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cyclopropyl = checker.check_ring("cyclopropane", mol_smiles)

            # If this is the target molecule (depth 0), check if it has cyclopropyl
            if depth == 0:
                target_has_cyclopropyl = has_cyclopropyl
                print(f"Target molecule: {mol_smiles}, Has cyclopropyl: {has_cyclopropyl}")

                # If target has cyclopropyl, add it to the set of molecules that should have it
                if has_cyclopropyl:
                    molecules_with_cyclopropyl.add(mol_smiles)
                else:
                    # If target doesn't have cyclopropyl, no need to check further
                    return

            # For non-target molecules, check if they should have cyclopropyl
            elif mol_smiles in molecules_with_cyclopropyl and not has_cyclopropyl:
                print(f"Molecule at depth {depth}: {mol_smiles}, Missing cyclopropyl")
                cyclopropyl_maintained = False

            # If this is a starting material, we don't need to check its children
            if node.get("in_stock", False):
                return

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, determine which reactants should have cyclopropyl
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # If the product should have cyclopropyl
            if product in molecules_with_cyclopropyl:
                # Check if this is a cyclopropyl-forming reaction
                product_has_cyclopropyl = checker.check_ring("cyclopropane", product)

                # Check if any reactant has cyclopropyl
                reactant_with_cyclopropyl = None
                for reactant in reactants:
                    if checker.check_ring("cyclopropane", reactant):
                        reactant_with_cyclopropyl = reactant
                        molecules_with_cyclopropyl.add(reactant)
                        break

                # If no reactant has cyclopropyl but product does, cyclopropyl was created in this reaction
                if not reactant_with_cyclopropyl and product_has_cyclopropyl:
                    print(f"Cyclopropyl created in reaction: {rsmi}")
                    cyclopropyl_creation_reactions.add(rsmi)
                    # Since cyclopropyl was created here, previous molecules don't need to have it
                    # We don't add any reactants to molecules_with_cyclopropyl
                elif reactant_with_cyclopropyl:
                    # If a reactant has cyclopropyl, it should be maintained
                    print(f"Cyclopropyl maintained from reactant: {reactant_with_cyclopropyl}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Only return True if target has cyclopropyl AND cyclopropyl is maintained
    result = target_has_cyclopropyl and cyclopropyl_maintained
    print(
        f"Final result: {result} (Target has cyclopropyl: {target_has_cyclopropyl}, Cyclopropyl maintained: {cyclopropyl_maintained})"
    )
    return result
