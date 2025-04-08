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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    Detects a synthetic strategy involving nucleophilic aromatic substitution to introduce an amine.
    """
    has_nucleophilic_aromatic_substitution = False

    def dfs_traverse(node):
        nonlocal has_nucleophilic_aromatic_substitution

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a nucleophilic aromatic substitution reaction
            if (
                checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
            ):
                print(f"Found nucleophilic aromatic substitution reaction: {rsmi}")
                has_nucleophilic_aromatic_substitution = True
            else:
                # Fallback check for nucleophilic aromatic substitution pattern
                has_aromatic_halide = False
                has_amine = False

                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aromatic_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        has_amine = True
                        print(f"Found amine in reactant: {reactant}")

                # Check if the product has a new C-N bond on an aromatic ring
                if has_aromatic_halide and has_amine:
                    p_mol = Chem.MolFromSmiles(product)
                    if p_mol and (
                        checker.check_fg("Aniline", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    ):
                        print(
                            f"Found C-N bond in product where halide was likely displaced: {product}"
                        )
                        has_nucleophilic_aromatic_substitution = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nucleophilic_aromatic_substitution
