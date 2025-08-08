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
    This function detects if the synthesis incorporates and maintains
    fluorinated aromatic groups throughout the synthesis.
    """
    # Track if the final product has fluorinated aromatic groups
    target_has_fluorinated = False

    # Track if fluorinated groups are maintained throughout synthesis
    fluorinated_maintained = True

    def has_fluorinated_aromatic(mol_smiles):
        """Helper function to check if a molecule has fluorinated aromatic groups"""
        # Check for aromatic halide
        has_aromatic_halide = checker.check_fg("Aromatic halide", mol_smiles)

        # If there's an aromatic halide, verify it's specifically fluorine
        if has_aromatic_halide:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F":
                        # Check if the fluorine is attached to an aromatic atom
                        neighbors = atom.GetNeighbors()
                        if neighbors and neighbors[0].GetIsAromatic():
                            return True

        # Also check for trifluoro groups
        return checker.check_fg("Trifluoro group", mol_smiles)

    def dfs_traverse(node, depth=0):
        nonlocal target_has_fluorinated, fluorinated_maintained

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for fluorinated aromatic groups
            has_fluoro_aromatic = has_fluorinated_aromatic(mol_smiles)

            # If this is the target molecule (root node), record if it has fluorinated aromatic groups
            if depth == 0:
                print(f"Target molecule: {mol_smiles}")
                if has_fluoro_aromatic:
                    print(f"Target molecule has fluorinated aromatic groups")
                    target_has_fluorinated = True
                else:
                    print(f"Target molecule does not have fluorinated aromatic groups")

            # For reactant molecules, check if they have fluorinated aromatic groups
            if "children" in node and node["children"]:
                for child in node["children"]:
                    if (
                        child["type"] == "reaction"
                        and "metadata" in child
                        and "rsmi" in child["metadata"]
                    ):
                        rxn_smiles = child["metadata"]["rsmi"]
                        product = rxn_smiles.split(">")[-1]
                        reactants = rxn_smiles.split(">")[0].split(".")

                        # Check if product has fluorinated aromatic groups
                        product_has_fluoro = has_fluorinated_aromatic(product)

                        # Check if any reactant has fluorinated aromatic groups
                        reactants_have_fluoro = any(
                            has_fluorinated_aromatic(reactant) for reactant in reactants
                        )

                        # In retrosynthesis, if product has fluorinated groups but reactants don't,
                        # this means fluorine was added in the forward direction
                        if product_has_fluoro and not reactants_have_fluoro:
                            print(
                                f"Fluorinated aromatic group added in forward reaction: {rxn_smiles}"
                            )
                            # This is fine for our strategy

                        # In retrosynthesis, if reactants have fluorinated groups but product doesn't,
                        # this means fluorine was removed in the forward direction
                        if reactants_have_fluoro and not product_has_fluoro:
                            print(
                                f"Fluorinated aromatic group removed in forward reaction: {rxn_smiles}"
                            )
                            fluorinated_maintained = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is successful if the target has fluorinated aromatic groups
    # and these groups are maintained throughout the synthesis
    result = target_has_fluorinated and fluorinated_maintained
    print(f"Strategy result: {result}")
    return result
