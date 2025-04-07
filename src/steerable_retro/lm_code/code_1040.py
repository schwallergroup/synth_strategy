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
    This function detects if the synthetic route involves construction of a purine scaffold.
    It looks for reactions where a purine ring system is formed from non-purine precursors.
    """
    purine_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal purine_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}: Checking reaction: {rsmi}")

                # Check if product contains purine-like structure
                from rdkit import Chem
                from rdkit.Chem import rdFMCS

                # First try the standard check
                product_has_purine = checker.check_ring("purine", product_smiles)

                # If standard check fails, try more detailed analysis
                if not product_has_purine:
                    # Define purine pattern
                    purine_pattern = "c1nc2ncnc2n1"
                    purine_mol = Chem.MolFromSmiles(purine_pattern)

                    # Clean product SMILES by removing atom mapping
                    clean_product = product_smiles
                    for i in range(30):  # Assuming max 30 atom mappings
                        clean_product = clean_product.replace(f":{i}]", "]")

                    product_mol = Chem.MolFromSmiles(clean_product)

                    if product_mol:
                        # Try to find purine substructure
                        if product_mol.HasSubstructMatch(purine_mol):
                            product_has_purine = True
                        else:
                            # Check for fused bicyclic system with correct number of nitrogens
                            ring_info = product_mol.GetRingInfo()
                            ring_atoms = ring_info.AtomRings()

                            # Look for atoms that are in multiple rings
                            atoms_in_multiple_rings = set()
                            for atom_idx in range(product_mol.GetNumAtoms()):
                                if ring_info.NumAtomRings(atom_idx) >= 2:
                                    atoms_in_multiple_rings.add(atom_idx)

                            # Check if we have a fused ring system with at least 4 nitrogens
                            if len(atoms_in_multiple_rings) >= 2:
                                n_count = 0
                                for ring in ring_atoms:
                                    for atom_idx in ring:
                                        if product_mol.GetAtomWithIdx(atom_idx).GetSymbol() == "N":
                                            n_count += 1

                                if n_count >= 4:  # Purine has 4 nitrogens
                                    product_has_purine = True

                print(f"Product contains purine: {product_has_purine}")

                if product_has_purine:
                    # Check if any reactant has purine
                    reactant_list = reactants_smiles.split(".")
                    reactants_have_purine = False

                    for r in reactant_list:
                        if checker.check_ring("purine", r):
                            reactants_have_purine = True
                            print(f"Reactant contains purine: {r}")
                            break

                        # If standard check fails, try more detailed analysis
                        if not reactants_have_purine:
                            # Clean reactant SMILES
                            clean_reactant = r
                            for i in range(30):
                                clean_reactant = clean_reactant.replace(f":{i}]", "]")

                            reactant_mol = Chem.MolFromSmiles(clean_reactant)

                            if reactant_mol and reactant_mol.HasSubstructMatch(purine_mol):
                                reactants_have_purine = True
                                print(f"Reactant contains purine (substructure): {r}")
                                break

                    if not reactants_have_purine:
                        print(f"Purine scaffold construction detected at depth {depth}")
                        purine_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: purine_formed = {purine_formed}")
    return purine_formed
