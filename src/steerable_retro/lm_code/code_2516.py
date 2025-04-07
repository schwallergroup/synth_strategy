#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects synthesis routes where specific aromatic motifs (chlorobenzene and phenyl ether)
    are maintained throughout the synthesis.
    """
    # Track if motifs are present in main pathway molecules
    main_pathway_molecules = []
    target_molecule = None

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            # Skip small molecules (likely reagents/solvents)
            if mol.GetNumAtoms() < 6:
                print(f"Skipping small molecule: {smiles}")
                return

            # The first molecule in DFS is the target molecule
            if depth == 0:
                nonlocal target_molecule
                target_molecule = {
                    "smiles": smiles,
                    "has_chlorobenzene": checker.check_fg("Aromatic halide", smiles),
                    "has_phenyl_ether": checker.check_fg("Ether", smiles)
                    and checker.check_ring("benzene", smiles),
                }
                print(f"Target molecule: {smiles}")
                print(f"  Has chlorobenzene: {target_molecule['has_chlorobenzene']}")
                print(f"  Has phenyl ether: {target_molecule['has_phenyl_ether']}")
            else:
                # For other molecules, check if they're part of the main pathway
                # by verifying they have similar structure to the target
                if mol.GetNumHeavyAtoms() > 10:  # Main pathway molecules are typically larger
                    main_pathway_molecules.append(
                        {
                            "smiles": smiles,
                            "has_chlorobenzene": checker.check_fg("Aromatic halide", smiles),
                            "has_phenyl_ether": checker.check_fg("Ether", smiles)
                            and checker.check_ring("benzene", smiles),
                        }
                    )
                    print(f"Main pathway molecule: {smiles}")
                    print(f"  Has chlorobenzene: {main_pathway_molecules[-1]['has_chlorobenzene']}")
                    print(f"  Has phenyl ether: {main_pathway_molecules[-1]['has_phenyl_ether']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have enough molecules in the main pathway
    if len(main_pathway_molecules) < 2:
        print(f"Not enough main pathway molecules: {len(main_pathway_molecules)}")
        return False

    # Check if target molecule has both motifs
    if not (target_molecule["has_chlorobenzene"] and target_molecule["has_phenyl_ether"]):
        print("Target molecule doesn't have both motifs")
        return False

    # Check if all main pathway molecules have both motifs
    all_have_chlorobenzene = all(mol["has_chlorobenzene"] for mol in main_pathway_molecules)
    all_have_phenyl_ether = all(mol["has_phenyl_ether"] for mol in main_pathway_molecules)

    result = all_have_chlorobenzene and all_have_phenyl_ether

    print(f"Maintained aromatic motifs strategy detected: {result}")
    print(f"All main pathway molecules have chlorobenzene: {all_have_chlorobenzene}")
    print(f"All main pathway molecules have phenyl ether: {all_have_phenyl_ether}")
    print(f"Total main pathway molecules: {len(main_pathway_molecules)}")

    return result
