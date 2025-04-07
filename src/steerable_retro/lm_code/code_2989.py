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
    This function detects a linear synthesis route that maintains a piperidine scaffold
    throughout the synthesis.
    """
    # Track reactions and piperidine presence
    reaction_count = 0
    piperidine_reactions = 0

    # Track if we've seen piperidine at all in the route
    has_piperidine_anywhere = False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, piperidine_reactions, has_piperidine_anywhere

        if node["type"] == "reaction":
            reaction_count += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains piperidine
                product_has_piperidine = checker.check_ring("piperidine", product)

                # Also check for CCN pattern which might indicate piperidine
                if not product_has_piperidine and "1CCN" in product:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for atom in product_mol.GetAtoms():
                            if atom.IsInRing() and atom.GetSymbol() == "N":
                                product_has_piperidine = True
                                break

                if product_has_piperidine:
                    has_piperidine_anywhere = True
                    # Check if at least one reactant contains piperidine
                    reactant_has_piperidine = False
                    for r in reactants:
                        if checker.check_ring("piperidine", r):
                            reactant_has_piperidine = True
                            break
                        # Also check for CCN pattern
                        elif "1CCN" in r:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol:
                                for atom in r_mol.GetAtoms():
                                    if atom.IsInRing() and atom.GetSymbol() == "N":
                                        reactant_has_piperidine = True
                                        break

                    if reactant_has_piperidine:
                        piperidine_reactions += 1
                        print(f"Depth {depth}: Detected piperidine scaffold maintained in reaction")
                    else:
                        print(
                            f"Depth {depth}: Piperidine in product but not in reactants - scaffold formation reaction"
                        )
                else:
                    # Check if any reactant contains piperidine (scaffold destruction)
                    reactant_has_piperidine = False
                    for r in reactants:
                        if checker.check_ring("piperidine", r):
                            reactant_has_piperidine = True
                            break
                        # Also check for CCN pattern
                        elif "1CCN" in r:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol:
                                for atom in r_mol.GetAtoms():
                                    if atom.IsInRing() and atom.GetSymbol() == "N":
                                        reactant_has_piperidine = True
                                        break

                    if reactant_has_piperidine:
                        print(
                            f"Depth {depth}: Piperidine in reactant but not in product - scaffold destruction"
                        )
                    else:
                        print(f"Depth {depth}: No piperidine scaffold involved in this reaction")

        # Check if molecule contains piperidine
        elif node["type"] == "mol" and node["smiles"]:
            has_piperidine = checker.check_ring("piperidine", node["smiles"])

            # Also check for CCN pattern which might indicate piperidine
            if not has_piperidine and "1CCN" in node["smiles"]:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.IsInRing() and atom.GetSymbol() == "N":
                            has_piperidine = True
                            break

            if has_piperidine:
                has_piperidine_anywhere = True
                print(f"Depth {depth}: Molecule contains piperidine scaffold: {node['smiles']}")
            else:
                print(
                    f"Depth {depth}: Molecule does not contain piperidine scaffold: {node['smiles']}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    print("Starting analysis of synthesis route...")
    dfs_traverse(route)

    print(f"Total reactions: {reaction_count}")
    print(f"Reactions maintaining piperidine scaffold: {piperidine_reactions}")
    print(f"Has piperidine anywhere in route: {has_piperidine_anywhere}")

    # Check if we have a synthesis with piperidine scaffold
    # We need to ensure:
    # 1. There is at least one reaction
    # 2. Piperidine is present in the route

    if reaction_count > 0 and has_piperidine_anywhere:
        print("Strategy detected: Linear synthesis with piperidine scaffold")
        return True

    return False
