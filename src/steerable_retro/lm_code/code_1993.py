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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthesis route involves heterocycle formation,
    by checking for an increase in the number of rings containing heteroatoms.
    """
    heterocycle_formation_found = False

    # List of common heterocycles to check
    heterocycle_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "thiophene",
        "thiopyran",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Heterocycle formation reaction detected: {reaction_type} at depth {depth}"
                        )
                        heterocycle_formation_found = True
                        return

                # Count specific heterocycles in product
                product_heterocycles = {}
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles[heterocycle] = (
                            product_heterocycles.get(heterocycle, 0) + 1
                        )

                # Count specific heterocycles in reactants
                reactant_heterocycles = {}
                for reactant in reactants:
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles[heterocycle] = (
                                reactant_heterocycles.get(heterocycle, 0) + 1
                            )

                # Check if any new heterocycle type appears in the product
                for heterocycle in product_heterocycles:
                    if heterocycle not in reactant_heterocycles:
                        print(f"New heterocycle type {heterocycle} formed at depth {depth}")
                        heterocycle_formation_found = True
                        return

                # Check if the count of any heterocycle type increases
                for heterocycle in product_heterocycles:
                    if product_heterocycles[heterocycle] > reactant_heterocycles.get(
                        heterocycle, 0
                    ):
                        print(f"Increased count of heterocycle {heterocycle} at depth {depth}")
                        heterocycle_formation_found = True
                        return

                # Fallback to general heterocycle counting if specific types aren't detected
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Count heterocycles in product
                    product_heterocycle_count = 0
                    ring_info = product_mol.GetRingInfo()
                    for ring_atoms in ring_info.AtomRings():
                        is_heterocycle = False
                        for atom_idx in ring_atoms:
                            atom = product_mol.GetAtomWithIdx(atom_idx)
                            if atom.GetAtomicNum() not in [1, 6]:  # Not H or C
                                is_heterocycle = True
                                break
                        if is_heterocycle:
                            product_heterocycle_count += 1

                    # Count heterocycles in reactants
                    reactant_heterocycle_count = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            ring_info = reactant_mol.GetRingInfo()
                            for ring_atoms in ring_info.AtomRings():
                                is_heterocycle = False
                                for atom_idx in ring_atoms:
                                    atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                    if atom.GetAtomicNum() not in [1, 6]:  # Not H or C
                                        is_heterocycle = True
                                        break
                                if is_heterocycle:
                                    reactant_heterocycle_count += 1

                    if product_heterocycle_count > reactant_heterocycle_count:
                        print(f"General heterocycle formation detected at depth {depth}")
                        heterocycle_formation_found = True
                        return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formation_found
