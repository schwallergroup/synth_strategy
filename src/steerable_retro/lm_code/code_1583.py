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
    Detects aromatic nucleophilic substitution where a chlorine on an aromatic ring
    is replaced by an amine group.
    """
    # Initialize tracking variable
    has_aromatic_nucleophilic_substitution = False

    def dfs_traverse(node):
        nonlocal has_aromatic_nucleophilic_substitution

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for specific aromatic nucleophilic substitution reactions
                if (
                    checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                ):

                    # Verify that we have an aromatic halide in reactants and an amine in product
                    has_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    )
                    has_amine_product = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    if has_aromatic_halide and has_amine_product:
                        print(f"Detected aromatic nucleophilic substitution in reaction: {rsmi}")
                        has_aromatic_nucleophilic_substitution = True

                # Alternative check for nucleophilic aromatic substitution
                # This covers cases where the specific reaction type might not be recognized
                elif any(checker.check_fg("Aromatic halide", r) for r in reactants_smiles):
                    # Check if any reactant has aromatic halide and product has amine
                    # but no reactant has the same amine group
                    has_amine_reactants = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    has_amine_product = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    if has_amine_product and not has_amine_reactants:
                        print(
                            f"Detected potential aromatic nucleophilic substitution in reaction: {rsmi}"
                        )
                        has_aromatic_nucleophilic_substitution = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Aromatic nucleophilic substitution strategy detected: {has_aromatic_nucleophilic_substitution}"
    )
    return has_aromatic_nucleophilic_substitution
