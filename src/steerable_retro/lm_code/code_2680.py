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
    Detects if the synthesis involves a dichlorophenyl group throughout.
    This means the final product contains a dichlorophenyl group, and this group
    is preserved throughout the synthesis route.
    """
    # Track if the final product has a dichlorophenyl group
    final_product_has_dichlorophenyl = False
    # Track if all reactions preserve the dichlorophenyl group
    dichlorophenyl_preserved = True

    def has_dichlorophenyl(smiles):
        """Helper function to check if a molecule has a dichlorophenyl group"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check if the molecule contains a benzene ring
        if not checker.check_ring("benzene", smiles):
            return False

        # Count chlorines attached to aromatic carbons
        aromatic_chlorine_count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 17:  # Chlorine
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIsAromatic():
                        aromatic_chlorine_count += 1

        # Return True if there are at least 2 chlorines attached to aromatic carbons
        return aromatic_chlorine_count >= 2

    def dfs_traverse(node, depth=0, path=None):
        nonlocal final_product_has_dichlorophenyl, dichlorophenyl_preserved

        if path is None:
            path = []

        # Check molecule nodes
        if node["type"] == "mol" and node.get("smiles"):
            smiles = node["smiles"]
            has_group = has_dichlorophenyl(smiles)

            # If this is the final product (depth 0), check if it has a dichlorophenyl group
            if depth == 0:
                final_product_has_dichlorophenyl = has_group
                if has_group:
                    print(f"Final product has dichlorophenyl group: {smiles}")
                else:
                    print(f"Final product does NOT have dichlorophenyl group: {smiles}")

            # Add to path for tracking
            path.append((node["type"], smiles, has_group))

        # Check reaction nodes
        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has dichlorophenyl
            product_has_group = has_dichlorophenyl(product)

            # Check if any reactant has dichlorophenyl
            reactant_has_group = any(has_dichlorophenyl(r) for r in reactants)

            # If a reactant has dichlorophenyl but the product doesn't, the group wasn't preserved
            if reactant_has_group and not product_has_group:
                dichlorophenyl_preserved = False
                print(f"Dichlorophenyl group not preserved in reaction: {rsmi}")

            # Add to path for tracking
            path.append((node["type"], rsmi, product_has_group))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path.copy())

    # Start traversal
    dfs_traverse(route)

    # The strategy is valid if the final product has a dichlorophenyl group
    # and this group is preserved throughout the synthesis
    result = final_product_has_dichlorophenyl and dichlorophenyl_preserved
    print(f"Strategy result: {result}")
    return result
