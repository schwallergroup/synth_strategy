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
    This function detects if the synthetic route follows a linear strategy with
    a late-stage amide bond disconnection.
    """
    late_stage_disconnection = False
    is_linear = True
    total_nodes = 0
    branching_nodes = 0

    # First pass to check if the route is mostly linear
    def check_linearity(node):
        nonlocal total_nodes, branching_nodes

        total_nodes += 1
        if len(node.get("children", [])) > 1:
            branching_nodes += 1

        for child in node.get("children", []):
            check_linearity(child)

    # Second pass to find late-stage amide disconnection
    def dfs_traverse(node, depth=0):
        nonlocal late_stage_disconnection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide bond disconnection at low depth (late stage)
            if depth <= 3:  # Late stage (depth 0, 1, 2, or 3)
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction using reaction checkers
                is_amide_coupling = any(
                    [
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        ),
                        checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi),
                        checker.check_reaction("Schotten-Baumann to ester", rsmi),
                        checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        ),
                        checker.check_reaction("Ester with primary amine to amide", rsmi),
                        checker.check_reaction("Ester with secondary amine to amide", rsmi),
                        checker.check_reaction("Ester with ammonia to amide", rsmi),
                        checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi),
                        checker.check_reaction("Acyl chloride with ammonia to amide", rsmi),
                        checker.check_reaction("Carboxylic acid to amide conversion", rsmi),
                    ]
                )

                # Check for azide reduction to amine (which can be used in amide formation)
                is_azide_reduction = any(
                    [
                        checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi),
                        checker.check_reaction("Reduction of nitro groups to amines", rsmi),
                    ]
                )

                # Check if product contains an amide group
                has_amide_product = any(
                    [
                        checker.check_fg("Primary amide", product),
                        checker.check_fg("Secondary amide", product),
                        checker.check_fg("Tertiary amide", product),
                    ]
                )

                if is_amide_coupling:
                    print(
                        f"Found late-stage amide bond disconnection via reaction check at depth {depth}"
                    )
                    late_stage_disconnection = True
                elif is_azide_reduction and has_amide_product:
                    print(
                        f"Found late-stage azide reduction to amine (for amide formation) at depth {depth}"
                    )
                    late_stage_disconnection = True
                # Alternatively, check for amide formation by examining reactants and product
                elif len(reactants) > 0 and has_amide_product:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol:
                        print(f"Product contains amide group: {product}")

                        # Check if we have appropriate reactants for amide formation
                        has_acid_derivative = False
                        has_amine_derivative = False
                        has_azide = False

                        for r in reactants:
                            # Check for acid derivatives
                            if any(
                                [
                                    checker.check_fg("Carboxylic acid", r),
                                    checker.check_fg("Acyl halide", r),
                                    checker.check_fg("Ester", r),
                                    checker.check_fg("Anhydride", r),
                                ]
                            ):
                                has_acid_derivative = True
                                print(f"Found acid derivative in reactants: {r}")

                            # Check for amine derivatives
                            if any(
                                [
                                    checker.check_fg("Primary amine", r),
                                    checker.check_fg("Secondary amine", r),
                                    checker.check_fg("Aniline", r),
                                    checker.check_fg("Ammonia", r),
                                ]
                            ):
                                has_amine_derivative = True
                                print(f"Found amine derivative in reactants: {r}")

                            # Check for azide (which can be reduced to amine)
                            if checker.check_fg("Azide", r):
                                has_azide = True
                                print(f"Found azide in reactants (potential amine source): {r}")

                        # Verify appropriate components are present
                        if (has_acid_derivative and has_amine_derivative) or has_azide:
                            print(
                                f"Found late-stage amide bond disconnection via reactant/product analysis at depth {depth}"
                            )
                            late_stage_disconnection = True

        # Process children nodes
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Check if the route is mostly linear
    check_linearity(route)
    linearity_ratio = 1.0 - (branching_nodes / max(1, total_nodes))
    is_linear = linearity_ratio >= 0.7  # At least 70% of nodes should be non-branching

    print(f"Total nodes: {total_nodes}, Branching nodes: {branching_nodes}")
    print(f"Linearity ratio: {linearity_ratio:.2f}")

    # Only proceed if the route is mostly linear
    if is_linear:
        print("Route is sufficiently linear, checking for late-stage amide disconnection")
        dfs_traverse(route)
    else:
        print(f"Route is not sufficiently linear (linearity ratio: {linearity_ratio:.2f})")

    result = is_linear and late_stage_disconnection
    print(
        f"Final result: {result} (is_linear: {is_linear}, late_stage_disconnection: {late_stage_disconnection})"
    )
    return result
