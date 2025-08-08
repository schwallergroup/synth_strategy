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
    Detects a strategy involving the formation of a bromine-containing heterocycle
    that is maintained throughout the synthesis.
    """
    # List of heterocycles to check
    heterocycle_types = [
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
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "furan",
        "thiophene",
    ]

    found_heterocycle_formation = False
    maintains_bromine = True
    heterocycle_with_bromine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation, maintains_bromine, heterocycle_with_bromine_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for bromine in product
                if (
                    not checker.check_fg("Aromatic halide", product_smiles)
                    and not checker.check_fg("Primary halide", product_smiles)
                    and not checker.check_fg("Secondary halide", product_smiles)
                    and not checker.check_fg("Tertiary halide", product_smiles)
                    and not checker.check_fg("Alkenyl halide", product_smiles)
                ):
                    if "Br" in product_smiles:  # Double-check with string search
                        print(f"Bromine present but not in recognized FG at depth {depth}")
                    else:
                        maintains_bromine = False
                        print(f"Bromine not maintained at depth {depth}")

                # Check for heterocycle formation (typically early in synthesis)
                if depth >= 3:  # Early in synthesis (higher depth)
                    # Count rings in reactants and product
                    reactant_has_heterocycle_with_br = False
                    for reactant in reactants_smiles:
                        for ring_type in heterocycle_types:
                            if checker.check_ring(ring_type, reactant) and "Br" in reactant:
                                reactant_has_heterocycle_with_br = True
                                break

                    product_has_heterocycle_with_br = False
                    for ring_type in heterocycle_types:
                        if checker.check_ring(ring_type, product_smiles) and "Br" in product_smiles:
                            product_has_heterocycle_with_br = True
                            break

                    # If product has heterocycle with Br but reactants don't, a heterocycle was formed
                    if product_has_heterocycle_with_br and not reactant_has_heterocycle_with_br:
                        found_heterocycle_formation = True
                        print(f"Found bromine-containing heterocycle formation at depth {depth}")

                # Check if any stage has a bromine-containing heterocycle
                if "Br" in product_smiles:
                    for ring_type in heterocycle_types:
                        if checker.check_ring(ring_type, product_smiles):
                            heterocycle_with_bromine_found = True
                            print(f"Bromine-containing heterocycle found at depth {depth}")
                            break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is valid if a bromine-containing heterocycle is formed and bromine is maintained
    return found_heterocycle_formation and maintains_bromine and heterocycle_with_bromine_found
