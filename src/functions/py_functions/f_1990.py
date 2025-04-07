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
    This function detects a convergent synthesis strategy where multiple heterocyclic fragments
    are coupled together, specifically looking for indole-pyrimidine coupling or key transformations
    on connected heterocycles.
    """
    has_heterocycle_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_heterocycle_coupling

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Extract depth information - try multiple formats
                depth = None
                # Try to extract from ID field
                id_str = node.get("metadata", {}).get("ID", "")
                depth_match = re.search(r"Depth:?\s*(\d+)", id_str)
                if depth_match:
                    depth = int(depth_match.group(1))
                else:
                    # Use traversal depth as fallback
                    depth = current_depth

                # Focus on early to mid-synthesis (depth 3-4)
                if depth in [3, 4]:
                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check for indole and pyrimidine in reactants and product
                    try:
                        product_has_indole = checker.check_ring("indole", product)
                        product_has_pyrimidine = checker.check_ring(
                            "pyrimidine", product
                        )

                        # Only proceed if product has both heterocycles
                        if product_has_indole and product_has_pyrimidine:
                            print(f"Product contains both indole and pyrimidine")

                            # Check which reactants contain which heterocycles
                            reactants_with_indole = []
                            reactants_with_pyrimidine = []

                            for i, reactant in enumerate(reactants):
                                if reactant.strip():  # Skip empty reactants
                                    try:
                                        if checker.check_ring("indole", reactant):
                                            reactants_with_indole.append(i)
                                        if checker.check_ring("pyrimidine", reactant):
                                            reactants_with_pyrimidine.append(i)
                                    except Exception as e:
                                        print(
                                            f"Error checking rings in reactant {i}: {e}"
                                        )

                            print(f"Reactants with indole: {reactants_with_indole}")
                            print(
                                f"Reactants with pyrimidine: {reactants_with_pyrimidine}"
                            )

                            # Check for convergent synthesis or key transformations
                            is_coupling = False
                            fg_transformation = False

                            # Check for coupling reactions if heterocycles come from different reactants
                            if (
                                reactants_with_indole
                                and reactants_with_pyrimidine
                                and not any(
                                    i in reactants_with_indole
                                    for i in reactants_with_pyrimidine
                                )
                            ):

                                coupling_reactions = [
                                    "Suzuki",
                                    "Buchwald-Hartwig",
                                    "N-arylation",
                                    "Stille",
                                    "Negishi",
                                    "Heck",
                                    "Sonogashira",
                                ]

                                for rxn_type in coupling_reactions:
                                    try:
                                        if checker.check_reaction(rxn_type, rsmi):
                                            print(
                                                f"Detected {rxn_type} coupling reaction"
                                            )
                                            is_coupling = True
                                            break
                                    except Exception as e:
                                        print(
                                            f"Error checking reaction {rxn_type}: {e}"
                                        )

                                # If specific coupling reaction not found, check for generic coupling
                                if not is_coupling:
                                    print("Detected heterocycle coupling (generic)")
                                    is_coupling = True

                            # Check for important functional group transformations on connected heterocycles
                            # Case 1: Both heterocycles already connected in reactant
                            if any(
                                i in reactants_with_indole
                                for i in reactants_with_pyrimidine
                            ):
                                print("Both heterocycles already connected in reactant")

                                # Check for sulfonamide to amine transformation
                                try:
                                    for i in range(len(reactants)):
                                        if (
                                            i in reactants_with_indole
                                            and i in reactants_with_pyrimidine
                                        ):
                                            if checker.check_fg(
                                                "Sulfonamide", reactants[i]
                                            ) and checker.check_fg(
                                                "Primary amine", product
                                            ):
                                                print(
                                                    "Detected sulfonamide to primary amine transformation"
                                                )
                                                fg_transformation = True
                                except Exception as e:
                                    print(
                                        f"Error checking sulfonamide transformation: {e}"
                                    )

                                # Check for specific sulfonamide deprotection reaction
                                try:
                                    if checker.check_reaction(
                                        "N-glutarimide deprotection", rsmi
                                    ) or checker.check_reaction(
                                        "Phthalimide deprotection", rsmi
                                    ):
                                        print("Detected amine deprotection reaction")
                                        fg_transformation = True
                                except Exception as e:
                                    print(f"Error checking deprotection reaction: {e}")

                                # Check for other relevant transformations
                                try:
                                    if checker.check_reaction(
                                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                        rsmi,
                                    ):
                                        print("Detected nitrogen acylation reaction")
                                        fg_transformation = True
                                except Exception as e:
                                    print(f"Error checking acylation reaction: {e}")

                                # Check for nucleophilic aromatic substitution
                                try:
                                    if checker.check_reaction(
                                        "nucl_sub_aromatic_ortho_nitro", rsmi
                                    ) or checker.check_reaction(
                                        "nucl_sub_aromatic_para_nitro", rsmi
                                    ):
                                        print(
                                            "Detected nucleophilic aromatic substitution"
                                        )
                                        fg_transformation = True
                                except Exception as e:
                                    print(
                                        f"Error checking nucleophilic substitution: {e}"
                                    )

                                # Check for reduction of nitro group
                                try:
                                    for i in range(len(reactants)):
                                        if (
                                            i in reactants_with_indole
                                            and i in reactants_with_pyrimidine
                                        ):
                                            if checker.check_fg(
                                                "Nitro group", reactants[i]
                                            ) and checker.check_fg(
                                                "Primary amine", product
                                            ):
                                                print(
                                                    "Detected nitro reduction to amine"
                                                )
                                                fg_transformation = True
                                except Exception as e:
                                    print(f"Error checking nitro reduction: {e}")

                                # General check for sulfonamide deprotection
                                try:
                                    if (
                                        any(
                                            "S(=O)" in reactant
                                            for reactant in reactants
                                        )
                                        and "NH2" in product
                                    ):
                                        print(
                                            "Detected potential sulfonamide deprotection (SMILES pattern)"
                                        )
                                        fg_transformation = True
                                except Exception as e:
                                    print(f"Error in general sulfonamide check: {e}")

                            # If either coupling or key transformation detected, mark as heterocycle coupling strategy
                            if is_coupling or fg_transformation:
                                print(
                                    f"Confirmed heterocycle coupling or key transformation at depth {depth}"
                                )
                                has_heterocycle_coupling = True
                    except Exception as e:
                        print(f"Error checking product rings: {e}")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return has_heterocycle_coupling
