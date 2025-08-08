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
    This function detects a synthetic strategy involving the synthesis of
    a heterocycle containing trifluoromethyl groups.
    """
    has_trifluoromethyl = False
    has_heterocycle = False
    trifluoromethyl_preserved = False

    # List of common heterocyclic rings to check
    heterocycle_rings = [
        "pyrazole",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl, has_heterocycle, trifluoromethyl_preserved

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for trifluoromethyl groups using the checker function
            if checker.check_fg("Trifluoro group", mol_smiles):
                has_trifluoromethyl = True
                print(f"Found trifluoromethyl group in molecule: {mol_smiles}")

            # Check for heterocycles using the checker function
            for ring in heterocycle_rings:
                if checker.check_ring(ring, mol_smiles):
                    has_heterocycle = True
                    print(f"Found heterocycle ({ring}) in molecule: {mol_smiles}")
                    break

            # Check if final product has trifluoromethyl
            if depth == 0 and checker.check_fg("Trifluoro group", mol_smiles):
                trifluoromethyl_preserved = True
                print(f"Trifluoromethyl group preserved in final product: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if required features are present
    strategy_present = has_trifluoromethyl and has_heterocycle and trifluoromethyl_preserved

    print(f"Trifluoromethyl-containing heterocycle synthesis strategy detected: {strategy_present}")
    print(f"  - Has trifluoromethyl group: {has_trifluoromethyl}")
    print(f"  - Has heterocycle: {has_heterocycle}")
    print(f"  - Trifluoromethyl preserved in final product: {trifluoromethyl_preserved}")

    return strategy_present
