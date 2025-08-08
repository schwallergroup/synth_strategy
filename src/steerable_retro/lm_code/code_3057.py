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
    """Check if the route involves Boc protection of an amine."""
    boc_protection_found = False

    def dfs(node, depth=0):
        nonlocal boc_protection_found

        if boc_protection_found:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            # Check for Boc protection reactions
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                boc_protection_found = True
                print(f"Found Boc protection reaction: {rsmi}")
                return

            # Check if product has Boc group
            product = rsmi.split(">")[-1]
            if "OC(=O)OC(C)(C)C" in product or "CC(C)(C)OC(=O)" in product:
                # Check if any reactant has a primary or secondary amine
                reactants = rsmi.split(">")[0].split(".")
                if any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    boc_protection_found = True
                    print(f"Found Boc protection by structure analysis: {rsmi}")
                    return

        # Also check if a molecule has a Boc-protected amine
        if node["type"] == "mol" and not node.get("in_stock", False):
            mol_smiles = node["smiles"]
            if "OC(=O)OC(C)(C)C" in mol_smiles or "CC(C)(C)OC(=O)" in mol_smiles:
                if "NC(=O)OC(C)(C)C" in mol_smiles or "CC(C)(C)OC(=O)N" in mol_smiles:
                    boc_protection_found = True
                    print(f"Found Boc-protected amine in molecule: {mol_smiles}")
                    return

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return boc_protection_found
