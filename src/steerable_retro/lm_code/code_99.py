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
    Detects a synthetic strategy involving N-demethylation of methylamino groups.
    """
    has_n_demethylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_demethylation

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1].split(".")

                # In retrosynthesis, the product is the starting material and reactants are the targets
                # So for N-demethylation, we're looking for N-methylation in the reverse direction

                # First check if this is a known N-methylation reaction
                methylation_reactions = [
                    "N-methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "DMS Amine methylation",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in methylation_reactions):
                    print(f"Found N-demethylation (reverse of methylation) at depth {depth}")
                    has_n_demethylation = True
                    return

                # Check for N-demethylation by examining functional group changes
                # In retrosynthesis: product (starting material) → reactants (targets)
                for p in product_smiles:
                    if len(p) <= 1:  # Skip very simple fragments
                        continue

                    for reactant in reactants_smiles:
                        if len(reactant) <= 1:  # Skip very simple fragments
                            continue

                        # Check for primary amine → secondary amine conversion (in retrosynthesis)
                        # This would be secondary → primary in forward synthesis (demethylation)
                        if (
                            checker.check_fg("Primary amine", p)
                            and checker.check_fg("Secondary amine", reactant)
                            and not checker.check_fg("Primary amine", reactant)
                        ):
                            print(
                                f"Found primary to secondary amine conversion (demethylation in forward) at depth {depth}"
                            )
                            has_n_demethylation = True
                            return

                        # Check for secondary amine → tertiary amine conversion (in retrosynthesis)
                        # This would be tertiary → secondary in forward synthesis (demethylation)
                        if (
                            checker.check_fg("Secondary amine", p)
                            and checker.check_fg("Tertiary amine", reactant)
                            and not checker.check_fg("Secondary amine", reactant)
                        ):
                            print(
                                f"Found secondary to tertiary amine conversion (demethylation in forward) at depth {depth}"
                            )
                            has_n_demethylation = True
                            return

                        # Additional check for reactions that might not be captured by the above
                        # Look for reactions that could involve demethylation
                        if (
                            checker.check_fg("Tertiary amine", p)
                            or checker.check_fg("Secondary amine", p)
                        ) and (
                            (
                                checker.check_fg("Secondary amine", reactant)
                                and not checker.check_fg("Tertiary amine", reactant)
                            )
                            or (
                                checker.check_fg("Primary amine", reactant)
                                and not checker.check_fg("Secondary amine", reactant)
                            )
                        ):
                            # This is a potential demethylation - check if it's specifically removing a methyl
                            # We're looking at the reverse reaction in retrosynthesis
                            print(f"Found potential N-demethylation at depth {depth}")
                            has_n_demethylation = True
                            return

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"N-demethylation detected: {has_n_demethylation}")

    return has_n_demethylation
