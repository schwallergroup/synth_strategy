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
    Detects if the synthesis route involves N-methylation of a secondary amine.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an N-methylation reaction
                methylation_reactions = [
                    "N-methylation",
                    "Methylation with MeI_secondary",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "DMS Amine methylation",
                    "Methylation",
                    "Parnes methylation",
                ]

                is_methylation = False
                for rxn_type in methylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_methylation = True
                        print(f"Found {rxn_type} reaction at depth {depth}")
                        break

                # If not identified as a specific methylation reaction, check for the characteristic transformation
                if not is_methylation:
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                    if (
                        has_secondary_amine
                        and has_tertiary_amine
                        and not any(checker.check_fg("Tertiary amine", r) for r in reactants)
                    ):
                        print(f"Found potential N-methylation reaction at depth {depth}")
                        is_methylation = True

                if is_methylation:
                    # Check if a secondary amine is being methylated
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                    if has_secondary_amine and has_tertiary_amine:
                        print(f"Confirmed N-methylation of secondary amine at depth {depth}")
                        for r in reactants:
                            if checker.check_fg("Secondary amine", r):
                                print(f"Reactant with secondary amine: {r}")
                        print(f"Product with tertiary amine: {product}")
                        result = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return result
