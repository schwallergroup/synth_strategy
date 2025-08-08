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
    Detects if the route contains a trifluoroethyl ether group (OCH2CF3) in the final product.
    Note: Despite the function name, this actually checks for trifluoroethyl ether, not trifluoromethyl ether (OCF3).
    """
    found_cf3_ether = False

    def dfs_traverse(node, current_depth=0):
        nonlocal found_cf3_ether

        if node["type"] == "mol" and "smiles" in node:
            if current_depth == 0:  # Final product
                smiles = node["smiles"]
                print(f"Checking final product: {smiles}")

                try:
                    # Check if the molecule contains a trifluoro group
                    has_trifluoro = checker.check_fg("Trifluoro group", smiles)
                    print(f"Has trifluoro group: {has_trifluoro}")

                    # Check if the molecule contains an ether
                    has_ether = checker.check_fg("Ether", smiles)
                    print(f"Has ether: {has_ether}")

                    if has_trifluoro and has_ether:
                        # Get the molecule to check if there's a trifluoroethyl ether group
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            # SMARTS pattern for trifluoroethyl ether: oxygen connected to CH2CF3
                            cf3_ether_pattern = Chem.MolFromSmarts("OCC(F)(F)F")
                            if mol.HasSubstructMatch(cf3_ether_pattern):
                                found_cf3_ether = True
                                print(
                                    f"Found trifluoroethyl ether group (OCH2CF3) in final product: {smiles}"
                                )
                            else:
                                print(
                                    "Trifluoro and ether groups exist but not as trifluoroethyl ether"
                                )
                except Exception as e:
                    print(f"Error checking for trifluoroethyl ether: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return found_cf3_ether
