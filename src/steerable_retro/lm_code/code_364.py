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
    This function detects a synthetic strategy involving the use of nitrile groups
    as key intermediates in heterocycle synthesis.
    """
    nitrile_intermediates = False
    heterocycle_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_intermediates, heterocycle_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for nitrile groups in reactants
            nitrile_in_reactants = False
            for reactant in reactants_smiles:
                if checker.check_fg("Nitrile", reactant):
                    nitrile_in_reactants = True
                    nitrile_intermediates = True
                    print(f"Nitrile found in reactant: {reactant}")
                    break

            # Check for heterocycle formation if nitrile is present in reactants
            if nitrile_in_reactants:
                # Check for specific nitrile-based heterocycle formation reactions
                if (
                    checker.check_reaction("tetrazole_terminal", rsmi)
                    or checker.check_reaction("tetrazole_connect_regioisomere_1", rsmi)
                    or checker.check_reaction("tetrazole_connect_regioisomere_2", rsmi)
                    or checker.check_reaction(
                        "Azide-nitrile click cycloaddition to tetrazole", rsmi
                    )
                    or checker.check_reaction("Azide-nitrile click cycloaddition to triazole", rsmi)
                    or checker.check_reaction("pyrazole", rsmi)
                    or checker.check_reaction("1,2,4-triazole_acetohydrazide", rsmi)
                    or checker.check_reaction("1,2,4-triazole_carboxylic-acid/ester", rsmi)
                    or checker.check_reaction("3-nitrile-pyridine", rsmi)
                    or checker.check_reaction("oxadiazole", rsmi)
                ):
                    heterocycle_formed = True
                    print(f"Nitrile-based heterocycle formation reaction detected: {rsmi}")

                # Also check if product contains heterocycles
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    heterocycle_rings = [
                        "tetrazole",
                        "triazole",
                        "pyrazole",
                        "oxadiazole",
                        "isoxazole",
                        "pyridine",
                    ]
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product_smiles):
                            print(f"Heterocycle {ring} found in product: {product_smiles}")
                            # Only mark as formed if we didn't already have this ring in the reactants
                            ring_in_reactants = False
                            for reactant in reactants_smiles:
                                if checker.check_ring(ring, reactant):
                                    ring_in_reactants = True
                                    break

                            if not ring_in_reactants:
                                heterocycle_formed = True
                                print(f"New heterocycle {ring} formed from nitrile")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if nitriles were used as intermediates and led to heterocycle formation
    result = nitrile_intermediates and heterocycle_formed
    print(f"Nitrile-based heterocycle synthesis strategy detected: {result}")
    return result
