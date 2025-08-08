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
    This function detects a synthetic strategy involving construction of a heterocyclic system,
    particularly focusing on indole-containing polycyclic systems.
    """
    heterocycle_construction_detected = False

    # List of heterocyclic rings to check
    heterocycles = [
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyridine",
        "pyrimidine",
    ]

    # List of heterocycle-forming reactions
    heterocycle_reactions = [
        "Fischer indole",
        "Paal-Knorr pyrrole",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_construction_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Product: {product}")
                print(f"Reactants: {reactants}")

                # Check for heterocycle-forming reactions
                for reaction_type in heterocycle_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Heterocycle-forming reaction detected: {reaction_type}")
                        heterocycle_construction_detected = True
                        return

                # Check if product contains a heterocyclic system
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check for heterocycles in product
                    product_has_heterocycle = False
                    product_heterocycles = []

                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product):
                            product_has_heterocycle = True
                            product_heterocycles.append(heterocycle)

                    if product_has_heterocycle:
                        print(f"Product contains heterocycles: {product_heterocycles}")

                        # Check if reactants don't have the same heterocyclic system
                        reactant_has_same_heterocycles = False

                        for reactant in reactants:
                            reactant_heterocycles = []
                            for heterocycle in product_heterocycles:
                                if checker.check_ring(heterocycle, reactant):
                                    reactant_heterocycles.append(heterocycle)

                            if set(reactant_heterocycles) == set(product_heterocycles):
                                reactant_has_same_heterocycles = True
                                print(f"Reactant {reactant} has the same heterocycles")
                                break

                        if not reactant_has_same_heterocycles:
                            print(f"Heterocycle construction detected at depth {depth}")
                            heterocycle_construction_detected = True

                            # Check ring count to identify polycyclic systems
                            ring_info = product_mol.GetRingInfo()
                            ring_count = ring_info.NumRings()
                            print(f"Product has {ring_count} rings")

                            # Additional check for polycyclic systems
                            if ring_count >= 2:
                                print(f"Polycyclic heterocycle construction detected")

        for child in node.get("children", []):
            if not heterocycle_construction_detected:  # Stop traversal if already detected
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_construction_detected
