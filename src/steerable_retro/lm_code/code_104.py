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
    This function detects N-alkylation in the late stage of synthesis (low depth).
    """
    late_stage_n_alkylation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_n_alkylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Consider depth 0-1 as late stage
            if current_depth <= 1:
                # Check if this is an N-alkylation reaction using the checker
                n_alkylation_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "N-methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                ]

                is_n_alkylation = any(
                    checker.check_reaction(rxn_type, rsmi) for rxn_type in n_alkylation_reactions
                )

                if is_n_alkylation:
                    print(f"Detected late-stage N-alkylation at depth {current_depth}")
                    late_stage_n_alkylation = True
                else:
                    # Fallback check using functional groups if reaction check fails
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for alkyl halide in reactants
                        has_alkyl_halide = any(
                            checker.check_fg(halide, r)
                            for r in reactants
                            for halide in ["Primary halide", "Secondary halide", "Tertiary halide"]
                        )

                        # Check for amine in reactants
                        has_amine = any(
                            checker.check_fg(amine, r)
                            for r in reactants
                            for amine in ["Primary amine", "Secondary amine"]
                        )

                        # Check for tertiary amine in product
                        has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                        if has_alkyl_halide and has_amine and has_tertiary_amine:
                            # Additional check: ensure the amine is actually alkylated
                            # This is a simplification - in a real scenario, you'd need to track atom mappings
                            print(
                                f"Detected late-stage N-alkylation (FG check) at depth {current_depth}"
                            )
                            late_stage_n_alkylation = True
                    except Exception as e:
                        print(f"Error in N-alkylation functional group check: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_n_alkylation
