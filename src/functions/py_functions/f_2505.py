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
    Detects if the synthetic route involves the formation of a nitrogen-containing
    heterocycle through cyclization.
    """
    found_heterocycle_formation = False

    # List of nitrogen-containing heterocycles to check
    n_heterocycles = [
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
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "benzotriazole",
    ]

    # List of heterocycle formation reactions to check
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
        "indole",
        "oxadiazole",
        "imidazole",
        "Pictet-Spengler",
    ]

    def count_rings(mol_smiles):
        """Count the number of rings in a molecule"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol:
            return mol.GetRingInfo().NumRings()
        return 0

    def contains_n_heterocycle(mol_smiles):
        """Check if molecule contains a nitrogen heterocycle"""
        for heterocycle in n_heterocycles:
            if checker.check_ring(heterocycle, mol_smiles):
                return True
        return False

    def is_heterocycle_formation_reaction(rxn_smiles):
        """Check if the reaction is a heterocycle formation reaction"""
        for reaction_type in heterocycle_formation_reactions:
            if checker.check_reaction(reaction_type, rxn_smiles):
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                reaction_smiles = node["metadata"]["rsmi"]

                # Check if this is a known heterocycle formation reaction
                if is_heterocycle_formation_reaction(reaction_smiles):
                    found_heterocycle_formation = True
                    print(f"Found heterocycle formation reaction: {reaction_smiles}")
                    return

                # If not a known reaction type, check for ring formation
                reactants_str, _, product_str = reaction_smiles.split(">")

                # Check if ring count increases in the product compared to reactants
                reactants = reactants_str.split(".")
                max_reactant_rings = max(
                    [count_rings(r) for r in reactants if r.strip()], default=0
                )
                product_rings = count_rings(product_str)

                # Check if nitrogen heterocycle is formed
                has_n_heterocycle_product = contains_n_heterocycle(product_str)

                # Check if any reactant already has the N-heterocycle
                has_n_heterocycle_reactant = any(
                    contains_n_heterocycle(r) for r in reactants if r.strip()
                )

                # Heterocycle formation: more rings in product, product has N-heterocycle,
                # and either reactants don't have N-heterocycle or product has more rings
                if (
                    product_rings > max_reactant_rings
                    and has_n_heterocycle_product
                    and (
                        not has_n_heterocycle_reactant
                        or product_rings > max_reactant_rings
                    )
                ):
                    found_heterocycle_formation = True
                    print(
                        f"Detected heterocycle formation through ring count analysis: {reaction_smiles}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if found_heterocycle_formation:
        print("Detected heterocycle formation strategy")
        return True
    return False
