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
    This function detects a synthesis route where a heterocycle (particularly imidazole)
    is formed in the late stages of the synthesis.
    """
    heterocycle_formation_depth = None
    max_depth = 0

    # List of heterocycles to check
    heterocycles = [
        "imidazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "pyrrole",
        "furan",
        "thiophene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # List of heterocycle formation reactions to check
    heterocycle_reactions = [
        "{imidazole}",
        "{triaryl-imidazole}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{pyrazole}",
        "{oxadiazole}",
        "{benzofuran}",
        "{benzothiophene}",
        "{indole}",
        "Paal-Knorr pyrrole synthesis",
        "Paal-Knorr pyrrole",
        "{Paal-Knorr pyrrole}",
        "Fischer indole",
        "{Fischer indole}",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is a heterocycle formation reaction
            is_heterocycle_formation = False

            # Method 1: Check using reaction types
            for reaction in heterocycle_reactions:
                if checker.check_reaction(reaction, rsmi):
                    print(f"Heterocycle formation reaction detected: {reaction}")
                    is_heterocycle_formation = True
                    break

            # Method 2: Check if product has heterocycle but reactants don't
            if not is_heterocycle_formation:
                product_has_heterocycle = False
                product_heterocycles = []

                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product_str):
                        product_has_heterocycle = True
                        product_heterocycles.append(heterocycle)
                        print(f"Product contains heterocycle: {heterocycle}")

                if product_has_heterocycle:
                    reactants = reactants_str.split(".")
                    reactant_heterocycles = set()

                    for reactant in reactants:
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                reactant_heterocycles.add(heterocycle)
                                print(f"Reactant contains heterocycle: {heterocycle}")

                    # Check if any heterocycle in product is not in reactants
                    new_heterocycles = [
                        h for h in product_heterocycles if h not in reactant_heterocycles
                    ]
                    if new_heterocycles:
                        is_heterocycle_formation = True
                        print(f"New heterocycles formed: {new_heterocycles}")

            if is_heterocycle_formation and heterocycle_formation_depth is None:
                heterocycle_formation_depth = depth
                print(f"Heterocycle formation detected at depth {depth}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if heterocycle formation occurs in the late stage (first half of synthesis)
    if heterocycle_formation_depth is not None and max_depth > 0:
        # Lower depth values correspond to later stages in synthesis
        is_late_stage = heterocycle_formation_depth <= (max_depth / 2)
        print(
            f"Heterocycle formation at depth {heterocycle_formation_depth}, max depth {max_depth}, is late stage: {is_late_stage}"
        )
        return is_late_stage

    return False
