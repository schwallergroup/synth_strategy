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
    This function detects if the synthesis uses silicon-containing protecting groups
    (like TMS, TBDMS, etc.) for alcohol protection.
    """
    silicon_pg_found = False

    def dfs_traverse(node, depth=0):
        nonlocal silicon_pg_found

        if node["type"] == "mol":
            # Check for silicon-containing protecting group in molecules
            if checker.check_fg("TMS ether protective group", node["smiles"]) or checker.check_fg(
                "Silyl protective group", node["smiles"]
            ):
                print(
                    f"Silicon protecting group detected in molecule at depth {depth}: {node['smiles']}"
                )
                silicon_pg_found = True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check for protection/deprotection reactions
            rxn_smiles = node["metadata"]["rsmi"]

            if checker.check_reaction("Alcohol protection with silyl ethers", rxn_smiles):
                print(f"Silicon protection reaction detected at depth {depth}: {rxn_smiles}")
                silicon_pg_found = True

            elif (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rxn_smiles)
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rxn_smiles
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rxn_smiles
                )
            ):
                print(f"Silicon deprotection reaction detected at depth {depth}: {rxn_smiles}")
                silicon_pg_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return silicon_pg_found
