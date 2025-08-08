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
    This function detects if the synthetic route employs a nucleophilic aromatic substitution strategy.
    It looks for reactions where a halogen on an aromatic ring is replaced by a nitrogen.
    """
    found_snar = False

    def dfs_traverse(node):
        nonlocal found_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a nucleophilic aromatic substitution reaction
            # First try to directly check the reaction type
            if (
                checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("N-arylation", rsmi)
                or checker.check_reaction("Buchwald-Hartwig", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
            ):
                print(f"Found nucleophilic aromatic substitution reaction: {rsmi}")
                found_snar = True
                return

            # If direct reaction check fails, look for characteristic patterns
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic halide in reactants
            has_aromatic_halide = False
            for reactant in reactants:
                try:
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aromatic_halide = True
                        break
                except Exception as e:
                    print(f"Error checking aromatic halide: {e}")
                    continue

            # Check for aniline or similar in product
            try:
                if has_aromatic_halide and (
                    checker.check_fg("Aniline", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                ):
                    # Additional check to ensure it's not just any reaction with these groups
                    # Look for nitrogen-containing heterocycles that might be formed
                    if (
                        checker.check_ring("pyridine", product)
                        or checker.check_ring("pyrimidine", product)
                        or checker.check_ring("pyrazine", product)
                        or checker.check_ring("indole", product)
                        or checker.check_ring("quinoline", product)
                        or checker.check_ring("isoquinoline", product)
                    ):
                        print(f"Found likely nucleophilic aromatic substitution reaction: {rsmi}")
                        found_snar = True
            except Exception as e:
                print(f"Error checking product: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_snar
