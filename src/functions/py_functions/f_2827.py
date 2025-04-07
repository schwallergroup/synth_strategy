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
    This function detects if the synthesis includes construction of a heterocyclic scaffold.
    It checks for the formation of various heterocyclic rings that weren't present in the reactants.
    """
    has_scaffold_construction = False

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "isothiazole",
        "pyrrole",
        "furan",
        "thiophene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
    ]

    # List of heterocycle-forming reaction types
    heterocycle_forming_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "oxadiazole",
        "imidazole",
    ]

    def dfs_traverse(node):
        nonlocal has_scaffold_construction

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known heterocycle-forming reaction
                for reaction_type in heterocycle_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected heterocycle-forming reaction: {reaction_type}")
                        has_scaffold_construction = True
                        return

                # Parse molecules
                reactant_mols = []
                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        reactant_mols.append(mol)

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and reactant_mols:
                    # Check for heterocyclic rings in product that aren't in reactants
                    for ring_name in heterocyclic_rings:
                        if checker.check_ring(ring_name, product):
                            # Check if any reactant already has this ring
                            reactants_have_ring = any(
                                checker.check_ring(ring_name, Chem.MolToSmiles(mol))
                                for mol in reactant_mols
                            )

                            if not reactants_have_ring:
                                print(
                                    f"Detected heterocyclic scaffold construction: {ring_name}"
                                )
                                has_scaffold_construction = True
                                return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_scaffold_construction
