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
    This function detects a convergent synthesis strategy where two complex fragments
    are combined in the final step.
    """
    convergent_synthesis_found = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check if we have at least 2 reactants
            if len(reactants) >= 2:
                # Check if reactants are complex enough
                reactant_mols = []
                complex_reactants = 0

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            reactant_mols.append(mol)
                            # Define "complex" as having more than 7 atoms and at least 1 ring
                            if (
                                mol.GetNumAtoms() > 7
                                and rdMolDescriptors.CalcNumRings(mol) > 0
                            ):
                                complex_reactants += 1
                    except Exception as e:
                        print(f"Error processing reactant: {e}")
                        continue

                # Check if product is formed
                try:
                    product_mol = Chem.MolFromSmiles(product_part)

                    # Convergent synthesis criteria:
                    # 1. At least 2 complex reactants
                    # 2. Check for common coupling reactions
                    is_coupling_reaction = (
                        checker.check_reaction(
                            "Suzuki coupling with boronic acids", rsmi
                        )
                        or checker.check_reaction(
                            "Suzuki coupling with boronic esters", rsmi
                        )
                        or checker.check_reaction("N-arylation", rsmi)
                        or checker.check_reaction("Buchwald-Hartwig", rsmi)
                        or checker.check_reaction("Sonogashira", rsmi)
                        or checker.check_reaction("Heck", rsmi)
                        or checker.check_reaction("Negishi coupling", rsmi)
                        or checker.check_reaction("Stille reaction", rsmi)
                    )

                    if complex_reactants >= 2 or (
                        complex_reactants >= 1 and is_coupling_reaction
                    ):
                        print(
                            f"Convergent synthesis detected with {complex_reactants} complex fragments"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Is coupling reaction: {is_coupling_reaction}")
                        convergent_synthesis_found = True

                except Exception as e:
                    print(f"Error processing product: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return convergent_synthesis_found
