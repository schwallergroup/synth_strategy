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
    This function detects the incorporation of multiple heterocyclic fragments
    in the synthesis.
    """
    # List of heterocycles to check for
    heterocycle_types = [
        "pyrimidine",
        "pyrazole",
        "pyridine",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # Track which heterocycles have been found
    found_heterocycles = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains any of our target heterocycles
            for heterocycle in heterocycle_types:
                if heterocycle not in found_heterocycles and checker.check_ring(
                    heterocycle, mol_smiles
                ):
                    found_heterocycles.add(heterocycle)
                    print(f"Detected {heterocycle} heterocycle at depth {depth}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                try:
                    # Extract reactants and product
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check reactants for heterocycles
                    for reactant_smiles in reactants_part.split("."):
                        for heterocycle in heterocycle_types:
                            if heterocycle not in found_heterocycles and checker.check_ring(
                                heterocycle, reactant_smiles
                            ):
                                found_heterocycles.add(heterocycle)
                                print(
                                    f"Detected {heterocycle} heterocycle in reactant at depth {depth}"
                                )

                    # Check product for heterocycles that might have been formed
                    for heterocycle in heterocycle_types:
                        if heterocycle not in found_heterocycles and checker.check_ring(
                            heterocycle, product_part
                        ):
                            # Verify this heterocycle wasn't in any of the reactants
                            reactant_has_heterocycle = any(
                                checker.check_ring(heterocycle, r)
                                for r in reactants_part.split(".")
                            )
                            if not reactant_has_heterocycle:
                                found_heterocycles.add(heterocycle)
                                print(
                                    f"Detected {heterocycle} heterocycle formation in product at depth {depth}"
                                )

                except Exception as e:
                    print(f"Error processing reaction SMILES at depth {depth}: {e}")

        # Process children (depth-first)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if at least 2 different heterocycles are incorporated
    print(f"Total unique heterocycles found: {len(found_heterocycles)}")
    print(f"Heterocycles found: {found_heterocycles}")
    return len(found_heterocycles) >= 2
