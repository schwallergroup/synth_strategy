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
    Detects a synthetic strategy involving ketone to amine conversion via oxime intermediate
    in the early stages of synthesis.
    """
    # Helper function to compare molecules ignoring atom mapping
    def compare_molecules_ignoring_mapping(smiles1, smiles2):
        try:
            # Remove atom mapping and convert to canonical SMILES
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol1 is None or mol2 is None:
                return False

            # Clear atom mapping
            for atom in mol1.GetAtoms():
                atom.SetAtomMapNum(0)
            for atom in mol2.GetAtoms():
                atom.SetAtomMapNum(0)

            # Get canonical SMILES without atom mapping
            canon_smiles1 = Chem.MolToSmiles(mol1, isomericSmiles=True)
            canon_smiles2 = Chem.MolToSmiles(mol2, isomericSmiles=True)

            print(f"Comparing canonical SMILES: {canon_smiles1} vs {canon_smiles2}")
            # Compare the canonical SMILES
            return canon_smiles1 == canon_smiles2
        except Exception as e:
            print(f"Error comparing molecules: {e}")
            return False

    # Track oxime formation and reduction with their depths
    oxime_formation = {"found": False, "depth": 0, "product_smiles": ""}
    oxime_reduction = {"found": False, "depth": 0, "reactant_smiles": ""}

    def dfs_traverse(node, depth=0):
        if node.get("type") == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxime formation
                if any(checker.check_fg("Ketone", r) for r in reactants) and checker.check_fg(
                    "Oxime", product
                ):
                    print(f"Found oxime formation at depth {depth}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")
                    oxime_formation["found"] = True
                    oxime_formation["depth"] = depth
                    oxime_formation["product_smiles"] = product

                # Check for oxime reduction
                if any(checker.check_fg("Oxime", r) for r in reactants) and checker.check_fg(
                    "Primary amine", product
                ):
                    print(f"Found oxime reduction at depth {depth}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")
                    oxime_reduction["found"] = True
                    oxime_reduction["depth"] = depth
                    for r in reactants:
                        if checker.check_fg("Oxime", r):
                            oxime_reduction["reactant_smiles"] = r
                            break

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both reactions were found and in early stages (depth >= 3)
    early_stage_threshold = 3
    sequence_found = (
        oxime_formation["found"]
        and oxime_reduction["found"]
        and oxime_formation["depth"] >= early_stage_threshold
        and oxime_reduction["depth"] >= early_stage_threshold
    )

    # Verify the oxime from formation is the same one used in reduction
    # Compare molecules ignoring atom mapping differences
    if sequence_found:
        oximes_match = compare_molecules_ignoring_mapping(
            oxime_formation["product_smiles"], oxime_reduction["reactant_smiles"]
        )

        if not oximes_match:
            print("Oxime formation and reduction are not connected - different molecules")
            # Try to print canonical SMILES for debugging
            try:
                form_mol = Chem.MolFromSmiles(oxime_formation["product_smiles"])
                red_mol = Chem.MolFromSmiles(oxime_reduction["reactant_smiles"])
                if form_mol and red_mol:
                    # Clear atom mapping
                    for atom in form_mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    for atom in red_mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    print(f"Formation oxime canonical: {Chem.MolToSmiles(form_mol)}")
                    print(f"Reduction oxime canonical: {Chem.MolToSmiles(red_mol)}")
            except Exception as e:
                print(f"Error generating canonical SMILES: {e}")

            sequence_found = False
        else:
            print("Oxime formation and reduction are connected - same molecule")

    print(f"Oxime formation: {oxime_formation}")
    print(f"Oxime reduction: {oxime_reduction}")
    print(f"Sequence found: {sequence_found}")

    return sequence_found
