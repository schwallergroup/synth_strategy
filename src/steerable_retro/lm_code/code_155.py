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
    This function detects a synthetic strategy involving thiazole and imidazole heterocycles.
    """
    has_thiazole = False
    has_imidazole = False
    has_heterocycle_formation = False

    # Additional nitrogen heterocycles that might be relevant
    nitrogen_heterocycles = ["pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole"]
    has_other_n_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_thiazole, has_imidazole, has_heterocycle_formation, has_other_n_heterocycle

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for thiazole
            if checker.check_ring("thiazole", mol_smiles):
                has_thiazole = True
                print(f"Detected thiazole heterocycle in molecule: {mol_smiles}")

            # Check for imidazole
            if checker.check_ring("imidazole", mol_smiles):
                has_imidazole = True
                print(f"Detected imidazole heterocycle in molecule: {mol_smiles}")

            # Check for other nitrogen heterocycles
            for ring in nitrogen_heterocycles:
                if checker.check_ring(ring, mol_smiles):
                    has_other_n_heterocycle = True
                    print(f"Detected {ring} heterocycle in molecule: {mol_smiles}")

        # Check for heterocycle formation reactions
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for thiazole formation
            if checker.check_reaction("thiazole", rxn_smiles):
                has_heterocycle_formation = True
                print(f"Detected thiazole formation reaction: {rxn_smiles}")

            # Check for imidazole formation
            if checker.check_reaction("imidazole", rxn_smiles):
                has_heterocycle_formation = True
                print(f"Detected imidazole formation reaction: {rxn_smiles}")

            # Check for benzimidazole formation (related to imidazole)
            if checker.check_reaction(
                "{benzimidazole_derivatives_carboxylic-acid/ester}", rxn_smiles
            ) or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rxn_smiles):
                has_heterocycle_formation = True
                print(f"Detected benzimidazole formation reaction: {rxn_smiles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Found thiazole: {has_thiazole}, Found imidazole: {has_imidazole}")
    print(f"Found other N-heterocycles: {has_other_n_heterocycle}")
    print(f"Found heterocycle formation reactions: {has_heterocycle_formation}")

    # Return True if thiazole OR imidazole is found, or if a heterocycle formation reaction is detected
    return has_thiazole or has_imidazole or has_heterocycle_formation
