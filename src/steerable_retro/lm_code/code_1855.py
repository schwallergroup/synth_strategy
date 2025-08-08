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
    Detects a strategy involving the synthesis of a compound with multiple
    different heteroaromatic systems (at least 3 different types).
    """
    # Track all heteroaromatic systems
    heteroaromatic_systems = {
        "pyridine": False,
        "pyrazole": False,
        "imidazole": False,
        "oxazole": False,
        "thiazole": False,
        "pyrimidine": False,
        "pyrazine": False,
        "pyridazine": False,
        "triazole": False,
        "tetrazole": False,
        "isoxazole": False,
        "isothiazole": False,
        "oxadiazole": False,
        "thiadiazole": False,
        "indole": False,
        "quinoline": False,
        "isoquinoline": False,
        "benzoxazole": False,
        "benzothiazole": False,
        "benzimidazole": False,
        "indazole": False,
        "benzotriazole": False,
        "furan": False,
        "thiophene": False,
    }

    def dfs_traverse(node):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for each heteroaromatic system
            for ring_name in heteroaromatic_systems.keys():
                if not heteroaromatic_systems[ring_name]:  # Only check if not already found
                    if checker.check_ring(ring_name, mol_smiles):
                        heteroaromatic_systems[ring_name] = True
                        print(f"Found {ring_name} in molecule: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Count how many different heteroaromatic systems we found
    heteroaromatic_count = sum(heteroaromatic_systems.values())
    result = heteroaromatic_count >= 3  # At least 3 different heteroaromatic systems

    print(f"Multi-heteroaromatic synthesis strategy detected: {result}")
    print(
        f"Found heteroaromatic systems: {[name for name, found in heteroaromatic_systems.items() if found]}"
    )
    print(f"Total heteroaromatic systems: {heteroaromatic_count}")

    return result
