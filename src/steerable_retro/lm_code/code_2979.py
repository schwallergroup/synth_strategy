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
    Detects a synthetic strategy involving amide bond disconnection.
    """
    has_amide_disconnection = False

    def dfs_traverse(node):
        nonlocal has_amide_disconnection

        if node["type"] == "reaction":
            metadata = node.get("metadata", {})
            rsmi = metadata.get("rsmi", "")
            if not rsmi:
                return

            # In retrosynthesis, products (in rsmi) are the starting materials
            # and reactants (in rsmi) are the target molecules
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation reactions directly
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
            ):
                print(f"Detected amide formation reaction: {rsmi}")
                has_amide_disconnection = True
                return

            # Alternative check: look for amide in product and carboxylic acid in reactants
            has_amide = (
                checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Secondary amide", product_smiles)
                or checker.check_fg("Tertiary amide", product_smiles)
            )

            has_carboxylic_acid = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
            )
            has_amine = any(checker.check_fg("Primary amine", r) for r in reactants_smiles) or any(
                checker.check_fg("Secondary amine", r) for r in reactants_smiles
            )

            if has_amide and (has_carboxylic_acid or has_amine):
                print(
                    f"Detected amide disconnection pattern: amide in product, acid/amine in reactants"
                )
                print(f"Reaction SMILES: {rsmi}")
                has_amide_disconnection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_amide_disconnection
