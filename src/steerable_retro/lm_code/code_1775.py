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
    Detects if the synthesis follows a linear strategy (no convergent steps).

    A linear synthesis builds the target molecule step by step, with each reaction
    adding complexity to a single growing intermediate. Convergent synthesis combines
    multiple complex intermediates in key steps.

    Returns:
        bool: True if the synthesis is linear, False if convergent steps are detected.
    """
    is_linear = True

    # Helper function to determine if a molecule is a simple reagent
    def is_simple_reagent(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return True  # Invalid SMILES treated as simple

            # Check for common reagents by functional groups
            if (
                checker.check_fg("Magnesium halide", smiles)
                or checker.check_fg("Alkyl lithium", smiles)
                or checker.check_fg("Aryl lithium", smiles)
                or checker.check_fg("Zinc halide", smiles)
                or checker.check_fg("Triflate", smiles)
                or checker.check_fg("Tosylate", smiles)
                or checker.check_fg("Mesylate", smiles)
            ):
                return True

            # Common solvents and bases
            common_solvents = [
                "CCO",
                "CCOC",
                "CC(=O)OC",
                "CS(=O)C",
                "CN(C)C=O",
                "c1ccccc1",
                "ClCCl",
                "CC(C)O",
                "O",
                "CO",
                "CC#N",
                "C1COCCO1",
                "CCOCC",
            ]
            if smiles in common_solvents:
                return True

            # Common catalysts and reagents often contain Pd, Pt, Rh, Ru, Cu
            if any(metal in smiles for metal in ["Pd", "Pt", "Rh", "Ru", "Cu", "Fe", "Ni", "Zn"]):
                return True

            # Simple molecules with few atoms
            if mol.GetNumAtoms() <= 8:
                return True

            # Simple molecules with low complexity
            ring_count = len(Chem.GetSSSR(mol))
            if mol.GetNumAtoms() <= 12 and ring_count == 0:
                return True

            # Common reagents
            if any(
                reagent in smiles
                for reagent in [
                    "B(O)O",
                    "P",
                    "S(=O)(=O)",
                    "C(=O)O",
                    "[Na+]",
                    "[K+]",
                    "[Li+]",
                    "Cl",
                    "Br",
                    "I",
                ]
            ):
                if mol.GetNumAtoms() <= 15:
                    return True

            return False
        except Exception as e:
            print(f"Error in is_simple_reagent: {e}")
            return True  # If error, assume it's simple

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:  # Skip if already found non-linear
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if this is a known reaction type that's inherently convergent
                # but typically used in linear strategies
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                ):
                    # These coupling reactions are often part of linear strategies
                    # despite having multiple reactants
                    pass
                else:
                    # Count significant reactants (not simple reagents)
                    significant_reactants = [r for r in reactants if not is_simple_reagent(r)]

                    if len(significant_reactants) > 1:
                        print(f"Detected convergent step at depth {depth}: {rsmi}")
                        print(f"Significant reactants: {significant_reactants}")
                        print(f"Number of significant reactants: {len(significant_reactants)}")
                        is_linear = False
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
