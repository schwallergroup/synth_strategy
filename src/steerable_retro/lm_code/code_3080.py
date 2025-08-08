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
    Detects a synthetic strategy involving sequential SNAr reactions
    on a halogenated heterocycle.
    """
    snar_reactions = []
    pyrimidine_molecules = []

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions, pyrimidine_molecules

        if node["type"] == "mol":
            # Check for pyrimidine core
            if checker.check_ring("pyrimidine", node["smiles"]):
                pyrimidine_molecules.append((node["smiles"], depth))
                print(f"Found pyrimidine at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr reaction
            is_snar = False

            # Check for nucleophilic aromatic substitution reactions
            if (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
            ):
                is_snar = True

            # If not found in reaction dictionary, check for characteristic patterns
            if not is_snar:
                # Check if any reactant has a pyrimidine ring
                has_pyrimidine_reactant = any(
                    checker.check_ring("pyrimidine", r) for r in reactants
                )

                # Check if any reactant has an aromatic halide
                has_aromatic_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                # Check for nucleophiles (amine or alcohol)
                has_nucleophile = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Phenol", r)
                    or checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    for r in reactants
                )

                # Check if product has a pyrimidine ring
                has_pyrimidine_product = checker.check_ring("pyrimidine", product)

                # If all conditions are met, it's likely an SNAr reaction
                if (
                    has_pyrimidine_reactant
                    and has_aromatic_halide
                    and has_nucleophile
                    and has_pyrimidine_product
                ):
                    is_snar = True

            if is_snar:
                snar_reactions.append(depth)
                print(f"Detected SNAr reaction at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present (at least 2 SNAr reactions and pyrimidine core)
    has_pyrimidine = len(pyrimidine_molecules) > 0
    strategy_present = has_pyrimidine and len(snar_reactions) >= 2

    print(f"Sequential SNAr strategy detection:")
    print(f"  Pyrimidine present: {has_pyrimidine}")
    print(f"  SNAr reactions: {len(snar_reactions)} at depths {snar_reactions}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
