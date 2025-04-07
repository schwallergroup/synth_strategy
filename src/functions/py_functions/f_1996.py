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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent),
    by checking if each reaction has only one complex reactant.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 6 atoms)
                complex_reactant_count = 0
                simple_reagent_count = 0

                for reactant in reactants:
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Consider a reactant complex if it has more than 6 atoms and is not a common reagent
                            is_common_reagent = (
                                checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                                or checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Tertiary amine", reactant)
                                or checker.check_fg("Carboxylic acid", reactant)
                            )

                            if reactant_mol.GetNumAtoms() > 6 and not is_common_reagent:
                                complex_reactant_count += 1
                                print(f"Complex reactant found: {reactant}")
                            else:
                                simple_reagent_count += 1
                                print(f"Simple reagent found: {reactant}")
                    except Exception as e:
                        print(f"Error processing reactant SMILES {reactant}: {e}")

                # A reaction is convergent if it has more than one complex reactant
                # We allow multiple simple reagents in a linear synthesis
                if complex_reactant_count > 1:
                    print(
                        f"Convergent step detected at depth {depth}: {complex_reactant_count} complex reactants"
                    )
                    is_linear = False
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
