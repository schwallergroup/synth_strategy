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
    Detects if the synthesis follows a linear strategy (no convergent steps with multiple complex fragments).
    """
    is_linear = True
    debug = False  # Set to True to enable debug prints

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is a known coupling reaction type
                is_coupling_reaction = any(
                    [
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi),
                        checker.check_reaction("Suzuki coupling with boronic esters", rsmi),
                        checker.check_reaction("Negishi coupling", rsmi),
                        checker.check_reaction("Stille reaction_aryl", rsmi),
                        checker.check_reaction("Heck terminal vinyl", rsmi),
                        checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi),
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        ),
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        ),
                    ]
                )

                # Count non-trivial reactants (more than 15 atoms)
                complex_reactants = 0
                for r in reactants:
                    if not r:
                        continue

                    mol = Chem.MolFromSmiles(r)
                    if not mol:
                        continue

                    # Skip common reagents by checking functional groups and size
                    is_common_reagent = False
                    is_small_molecule = mol.GetNumAtoms() <= 8

                    # Check for common protecting groups and reagents
                    if (
                        checker.check_fg("Boc", r)
                        or checker.check_fg("TMS ether protective group", r)
                        or checker.check_fg("Silyl protective group", r)
                        or checker.check_fg("Triflate", r)
                        or checker.check_fg("Tosylate", r)
                        or checker.check_fg("Mesylate", r)
                        or checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Aromatic halide", r)
                    ):
                        is_common_reagent = True

                    # Check for common coupling reagents
                    if (
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        or checker.check_fg("Magnesium halide", r)
                        or checker.check_fg("Zinc halide", r)
                        or checker.check_fg("Tin", r)
                        or checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Alkyne", r)
                    ):
                        is_common_reagent = True

                    # Count as complex only if it has >15 atoms and is not a common reagent or small molecule
                    if mol.GetNumAtoms() > 15 and not is_common_reagent and not is_small_molecule:
                        complex_reactants += 1
                        if debug:
                            print(f"Complex reactant found: {r} with {mol.GetNumAtoms()} atoms")

                # If any step has more than one complex reactant, or is a coupling reaction with multiple reactants,
                # it's not a linear synthesis
                if complex_reactants > 1 or (
                    is_coupling_reaction and len(reactants) > 1 and complex_reactants > 0
                ):
                    is_linear = False
                    if debug:
                        print(f"Found convergent step with multiple complex reactants: {rsmi}")
                        if is_coupling_reaction:
                            print(f"This is a coupling reaction: {rsmi}")
            except Exception as e:
                if debug:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
