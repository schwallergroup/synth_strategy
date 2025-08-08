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
    Detects if the synthesis route involves protection of a carboxylic acid
    in the late stage (low depth) of the synthesis.
    """
    protection_at_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal protection_at_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for protection reaction (forward direction)
                has_acid_in_reactants = any(
                    checker.check_fg("Carboxylic acid", smi) for smi in reactants_smiles
                )
                has_ester_in_product = checker.check_fg("Ester", product_smiles)

                # Check for deprotection reaction (forward direction, which would be protection in retrosynthesis)
                has_ester_in_reactants = any(
                    checker.check_fg("Ester", smi) for smi in reactants_smiles
                )
                has_acid_in_product = checker.check_fg("Carboxylic acid", product_smiles)

                # Check if this is a protection or deprotection reaction
                is_protection_reaction = (
                    checker.check_reaction("Protection of carboxylic acid", rsmi)
                    or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                )

                is_deprotection_reaction = (
                    checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or (
                        has_ester_in_reactants and has_acid_in_product
                    )  # Fallback check for deprotection
                )

                print(f"  Has acid in reactants: {has_acid_in_reactants}")
                print(f"  Has ester in product: {has_ester_in_product}")
                print(f"  Has ester in reactants: {has_ester_in_reactants}")
                print(f"  Has acid in product: {has_acid_in_product}")
                print(f"  Is protection reaction: {is_protection_reaction}")
                print(f"  Is deprotection reaction: {is_deprotection_reaction}")

                # Forward direction protection
                if has_acid_in_reactants and has_ester_in_product and is_protection_reaction:
                    # Verify that the carboxylic acid is being protected (not just present)
                    acid_count_reactants = sum(
                        len(checker.get_fg_atom_indices("Carboxylic acid", smi))
                        for smi in reactants_smiles
                    )
                    acid_count_product = len(
                        checker.get_fg_atom_indices("Carboxylic acid", product_smiles)
                    )

                    print(f"  Acid count in reactants: {acid_count_reactants}")
                    print(f"  Acid count in product: {acid_count_product}")

                    if acid_count_reactants > acid_count_product:
                        protection_at_depth = depth
                        print(f"Detected carboxylic acid protection at depth {depth}")

                # Reverse direction protection (deprotection in forward direction)
                elif has_ester_in_reactants and has_acid_in_product and is_deprotection_reaction:
                    # In retrosynthesis, this would be a protection
                    ester_count_reactants = sum(
                        len(checker.get_fg_atom_indices("Ester", smi)) for smi in reactants_smiles
                    )
                    ester_count_product = len(checker.get_fg_atom_indices("Ester", product_smiles))
                    acid_count_reactants = sum(
                        len(checker.get_fg_atom_indices("Carboxylic acid", smi))
                        for smi in reactants_smiles
                    )
                    acid_count_product = len(
                        checker.get_fg_atom_indices("Carboxylic acid", product_smiles)
                    )

                    print(f"  Ester count in reactants: {ester_count_reactants}")
                    print(f"  Ester count in product: {ester_count_product}")
                    print(f"  Acid count in reactants: {acid_count_reactants}")
                    print(f"  Acid count in product: {acid_count_product}")

                    if ester_count_reactants > 0 and acid_count_product > 0:
                        # This is a deprotection in forward direction, which means protection in retrosynthesis
                        protection_at_depth = depth
                        print(
                            f"Detected carboxylic acid protection (from deprotection) at depth {depth}"
                        )

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it happens at depth 0 or 1
    is_late_stage = protection_at_depth is not None and protection_at_depth <= 1

    if is_late_stage:
        print(f"Confirmed late-stage carboxylic acid protection at depth {protection_at_depth}")
    else:
        if protection_at_depth is not None:
            print(
                f"Carboxylic acid protection found, but not in late stage (depth {protection_at_depth})"
            )
        else:
            print("No carboxylic acid protection detected in the synthesis route")

    return is_late_stage
