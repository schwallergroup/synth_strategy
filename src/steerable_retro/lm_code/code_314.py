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
    This function detects a mid-stage N-alkylation reaction.
    Mid-stage is defined as occurring at depth 1 in the synthesis tree.
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction" and depth == 1:  # Mid stage
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Use checker functions to identify N-alkylation reactions
                is_primary_alkylation = checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rsmi
                )
                is_secondary_alkylation = checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                is_epoxide_ring_opening = checker.check_reaction(
                    "Ring opening of epoxide with amine", rsmi
                )

                # Check for other N-alkylation reaction types
                is_methylation_primary = checker.check_reaction(
                    "Methylation with MeI_primary", rsmi
                )
                is_methylation_secondary = checker.check_reaction(
                    "Methylation with MeI_secondary", rsmi
                )
                is_n_methylation = checker.check_reaction("N-methylation", rsmi)
                is_eschweiler_clarke_primary = checker.check_reaction(
                    "Eschweiler-Clarke Primary Amine Methylation", rsmi
                )
                is_eschweiler_clarke_secondary = checker.check_reaction(
                    "Eschweiler-Clarke Secondary Amine Methylation", rsmi
                )
                is_reductive_methylation = checker.check_reaction(
                    "Reductive methylation of primary amine with formaldehyde", rsmi
                )

                print(f"Primary alkylation: {is_primary_alkylation}")
                print(f"Secondary alkylation: {is_secondary_alkylation}")
                print(f"Epoxide ring opening: {is_epoxide_ring_opening}")
                print(
                    f"Methylation reactions: {is_methylation_primary}, {is_methylation_secondary}, {is_n_methylation}"
                )
                print(
                    f"Eschweiler-Clarke: {is_eschweiler_clarke_primary}, {is_eschweiler_clarke_secondary}"
                )
                print(f"Reductive methylation: {is_reductive_methylation}")

                if (
                    is_primary_alkylation
                    or is_secondary_alkylation
                    or is_methylation_primary
                    or is_methylation_secondary
                    or is_n_methylation
                    or is_eschweiler_clarke_primary
                    or is_eschweiler_clarke_secondary
                    or is_reductive_methylation
                    or is_epoxide_ring_opening
                ):
                    print(f"Detected N-alkylation reaction at depth {depth}")
                    n_alkylation_detected = True
                else:
                    # Fallback check using reactants and products
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for amine in reactants
                        amine_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                amine_in_reactants = True
                                print(f"Found amine in reactant: {reactant}")
                                break

                        # Check for tertiary amine in product
                        tertiary_amine_in_product = checker.check_fg("Tertiary amine", product)
                        if tertiary_amine_in_product:
                            print(f"Found tertiary amine in product: {product}")

                        # Check for alkyl halide in reactants
                        alkyl_halide_in_reactants = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary halide", reactant)
                                or checker.check_fg("Secondary halide", reactant)
                                or checker.check_fg("Tertiary halide", reactant)
                            ):
                                alkyl_halide_in_reactants = True
                                print(f"Found alkyl halide in reactant: {reactant}")
                                break

                        # Check for epoxide in reactants
                        epoxide_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Oxirane", reactant) or checker.check_ring(
                                "oxirane", reactant
                            ):
                                epoxide_in_reactants = True
                                print(f"Found epoxide in reactant: {reactant}")
                                break

                        if (
                            amine_in_reactants
                            and tertiary_amine_in_product
                            and (alkyl_halide_in_reactants or epoxide_in_reactants)
                        ):
                            print(
                                f"Detected N-alkylation through functional group analysis at depth {depth}"
                            )
                            n_alkylation_detected = True
                    except Exception as e:
                        print(f"Error in fallback check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return n_alkylation_detected
