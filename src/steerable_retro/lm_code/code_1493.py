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
    This function detects if the synthesis follows a linear (non-convergent) approach.

    A linear synthesis is characterized by each reaction step having only one major
    reactant (the previous intermediate) plus reagents, as opposed to convergent
    synthesis where multiple complex intermediates are combined.
    """
    is_linear = True

    # Define common reagent patterns
    common_reagents = [
        "Grignard",
        "Wittig",
        "Suzuki",
        "Buchwald-Hartwig",
        "Heck",
        "Sonogashira",
        "Stille",
        "Negishi",
        "Kumada",
        "Chan-Lam",
        "Mitsunobu",
        "Acylation",
        "Alkylation",
        "Reduction",
        "Oxidation",
    ]

    def is_likely_reagent(smiles):
        """Determine if a molecule is likely a reagent rather than a major reactant."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Small molecules are often reagents
            if mol.GetNumAtoms() < 8:
                return True

            # Check for common reagent functional groups
            for fg in [
                "Triflate",
                "Tosylate",
                "Mesylate",
                "Magnesium halide",
                "Acyl halide",
                "Boronic acid",
                "Boronic ester",
            ]:
                if checker.check_fg(fg, smiles):
                    return True

            # Common solvents and bases
            common_reagent_smiles = [
                "O",
                "CCO",
                "CCCO",
                "CC(C)O",
                "CN(C)C=O",
                "CC(=O)O",
                "CS(=O)(=O)O",
                "CC#N",
                "ClCCl",
                "C1CCOC1",
                "C1COCCO1",
            ]
            if smiles in common_reagent_smiles:
                return True

            return False
        except:
            # If analysis fails, default to not a reagent
            return False

    def count_major_reactants(reactants_list):
        """Count the number of major reactants (non-reagents) in a reaction."""
        major_count = 0
        for r in reactants_list:
            if not is_likely_reagent(r):
                major_count += 1
        return major_count

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:  # Early return if we already found it's not linear
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Check reaction type for special handling
            reaction_type = None
            for rxn in common_reagents:
                if (
                    "metadata" in node
                    and "ID" in node["metadata"]
                    and rxn in node["metadata"]["ID"]
                ):
                    reaction_type = rxn
                    break

            # Count major reactants
            major_reactant_count = count_major_reactants(reactants)

            # Special handling for coupling reactions which are inherently convergent
            # but might be part of a linear strategy
            if reaction_type in [
                "Suzuki",
                "Buchwald-Hartwig",
                "Heck",
                "Sonogashira",
                "Stille",
                "Negishi",
                "Kumada",
                "Chan-Lam",
            ]:
                # These are typically considered linear despite having two major components
                pass
            elif major_reactant_count > 1:
                print(
                    f"Found convergent step with {major_reactant_count} major reactants at depth {depth}"
                )
                print(f"Reaction SMILES: {rsmi}")
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Synthesis follows a linear (non-convergent) approach")
    else:
        print("Synthesis follows a convergent approach")

    return is_linear
