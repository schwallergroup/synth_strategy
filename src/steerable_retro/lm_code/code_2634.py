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
    Detects if the final product contains 3 or more aromatic rings.
    """
    aromatic_ring_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_ring_count

        if node["type"] == "mol" and depth == 0:  # Final product
            print(f"Analyzing final product: {node['smiles']}")
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                print("Failed to parse molecule SMILES")
                return

            # Method 1: Use RingInfo to find aromatic rings
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
            aromatic_rings = []

            for ring in atom_rings:
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    aromatic_rings.append(ring)
                    print(f"Found aromatic ring: {ring}")

            aromatic_ring_count = len(aromatic_rings)

            # Method 2: Use ring pattern detection as backup
            if aromatic_ring_count < 3:
                print("Using ring pattern detection as backup")
                aromatic_systems = [
                    "benzene",
                    "pyridine",
                    "pyrrole",
                    "furan",
                    "thiophene",
                    "imidazole",
                    "pyrazole",
                    "oxazole",
                    "thiazole",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                    "naphthalene",
                    "anthracene",
                ]

                # Use a set to avoid double-counting rings
                detected_rings = set()
                for ring_name in aromatic_systems:
                    if checker.check_ring(ring_name, node["smiles"]):
                        ring_indices = checker.get_ring_atom_indices(ring_name, node["smiles"])
                        print(f"Detected {len(ring_indices)} instances of {ring_name}")
                        for ring_atoms in ring_indices:
                            # Add tuple of sorted atom indices to ensure uniqueness
                            detected_rings.add(tuple(sorted(ring_atoms[0])))

                # Use the higher count between the two methods
                backup_count = len(detected_rings)
                print(f"Backup method found {backup_count} unique aromatic rings")
                aromatic_ring_count = max(aromatic_ring_count, backup_count)

            print(f"Total aromatic rings detected: {aromatic_ring_count}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = aromatic_ring_count >= 3
    print(f"3+ aromatic rings detected: {result}")
    return result
