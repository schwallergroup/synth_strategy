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
    Detects a synthetic route with late-stage coupling of two complex fragments.
    """
    has_late_stage_coupling = False

    # List of common coupling reactions
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Negishi coupling",
        "Heck terminal vinyl",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Ullmann condensation",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "decarboxylative_coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_coupling

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check at depths 0 and 1 (final and penultimate steps)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least 2 reactants (fragments)
            if len(reactants_smiles) >= 2:
                reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if both fragments are complex (have at least 10 atoms)
                complex_fragments = [mol for mol in reactants if mol and mol.GetNumAtoms() >= 10]

                if len(complex_fragments) >= 2:
                    # Check if this is a known coupling reaction
                    is_coupling_reaction = False
                    for reaction_type in coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Found coupling reaction: {reaction_type}")
                            is_coupling_reaction = True
                            break

                    # If it's a coupling reaction or we can verify fragments are joined
                    if is_coupling_reaction:
                        print(f"Found late-stage coupling of complex fragments at depth {depth}")
                        has_late_stage_coupling = True
                    else:
                        # Alternative check: verify that the product has more atoms than each individual reactant
                        # This suggests fragments were joined
                        if product and all(
                            product.GetNumAtoms() > mol.GetNumAtoms() for mol in reactants
                        ):
                            print(f"Found late-stage fragment joining at depth {depth}")
                            has_late_stage_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage fragment coupling detected: {has_late_stage_coupling}")
    return has_late_stage_coupling
