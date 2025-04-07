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
    This function detects if the synthetic route follows a linear synthesis strategy
    without convergent steps.

    A linear synthesis strategy is characterized by:
    1. Each reaction step has only one significant reactant (excluding reagents/solvents)
    2. The route doesn't branch into multiple parallel paths
    3. Certain reaction types (protection, deprotection, etc.) are considered linear
       even if they have multiple reactants
    """
    is_linear = True

    # List of reaction types that are inherently linear despite multiple reactants
    linear_reaction_types = [
        # Protection/deprotection reactions
        "Boc amine protection",
        "Boc amine deprotection",
        "Alcohol protection with silyl ethers",
        "Alcohol deprotection from silyl ethers",
        "Protection of carboxylic acid",
        "Deprotection of carboxylic acid",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "TMS deprotection from alkyne",
        "Tert-butyl deprotection of amine",
        # Common coupling reactions used in linear synthesis
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Heck terminal vinyl",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        # Other common linear synthesis reactions
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Esterification of Carboxylic Acids",
        "Williamson Ether Synthesis",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
    ]

    def is_common_reagent(smiles):
        """Check if a molecule is likely a common reagent rather than a significant reactant"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Common reagents typically have fewer heavy atoms
            if mol.GetNumHeavyAtoms() <= 5:
                return True

            # Check for common reagent functional groups
            if (
                checker.check_fg("Triflate", smiles)
                or checker.check_fg("Tosylate", smiles)
                or checker.check_fg("Mesylate", smiles)
            ):
                return True

            # Check for simple alcohols, acids, etc.
            if smiles in [
                "O",
                "OO",
                "CO",
                "CCO",
                "[OH-]",
                "[H+]",
                "O=C=O",
                "CC(=O)O",
                "C(=O)O",
                "CN",
            ]:
                return True

            return False
        except:
            return False

    def is_boronic_derivative(smiles):
        """Check if a molecule is a boronic acid or ester"""
        try:
            return checker.check_fg("Boronic acid", smiles) or checker.check_fg(
                "Boronic ester", smiles
            )
        except:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Check for branching structure (multiple children for a molecule node)
        if node["type"] == "mol" and len(node.get("children", [])) > 1:
            print(f"Branching detected at depth {depth} - molecule has multiple reaction paths")
            is_linear = False
            return

        if node["type"] == "reaction":
            # Check number of reactants in this reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if this is a known linear reaction type
                is_linear_reaction_type = False
                for rxn_type in linear_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Linear reaction type detected: {rxn_type}")
                        is_linear_reaction_type = True
                        break

                # Check for Suzuki-type coupling if not already identified
                if not is_linear_reaction_type:
                    # Count boronic derivatives - common in linear synthesis coupling reactions
                    boronic_count = sum(1 for r in reactants if is_boronic_derivative(r))
                    if boronic_count > 0 and len(reactants) <= 3:
                        print(
                            f"Potential coupling reaction with boronic derivative detected at depth {depth}"
                        )
                        is_linear_reaction_type = True

                if not is_linear_reaction_type:
                    # Filter out common reagents
                    significant_reactants = [r for r in reactants if not is_common_reagent(r)]

                    # If we have more than one significant reactant, it might be convergent
                    if len(significant_reactants) > 1:
                        # Additional check for reactions that might be linear despite multiple reactants
                        # For example, if one reactant is much larger than others, it might be the main scaffold
                        main_scaffold_found = False
                        for r in significant_reactants:
                            mol = Chem.MolFromSmiles(r)
                            if (
                                mol and mol.GetNumHeavyAtoms() > 15
                            ):  # Larger threshold for main scaffold
                                main_scaffold_found = True
                                break

                        if not main_scaffold_found and len(significant_reactants) > 1:
                            print(
                                f"Convergent step detected with multiple significant reactants at depth {depth}"
                            )
                            for r in significant_reactants:
                                print(f"Reactant SMILES: {r}")
                            is_linear = False
            else:
                print(f"Warning: Reaction node at depth {depth} missing metadata or rsmi")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
