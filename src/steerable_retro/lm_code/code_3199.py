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
    This function detects if the synthesis involves a late-stage coupling
    with a sulfonamide fragment.
    """
    found_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Late stage (up to 2 steps from final product)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if product contains sulfonamide
                sulfonamide_in_product = checker.check_fg("Sulfonamide", product)
                print(f"Sulfonamide in product: {sulfonamide_in_product}")

                if sulfonamide_in_product:
                    # Check if it's a known sulfonamide synthesis reaction
                    is_sulfonamide_rxn = checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ) or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )

                    # Check if any reactant contains sulfonyl halide (precursor to sulfonamide)
                    sulfonyl_halide_in_reactants = any(
                        checker.check_fg("Sulfonyl halide", reactant) for reactant in reactants
                    )

                    # Check if any reactant already contains sulfonamide (for coupling reactions)
                    sulfonamide_in_reactants = any(
                        checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                    )

                    # Check if any reactant contains amine (needed for sulfonamide formation)
                    amine_in_reactants = any(
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        for reactant in reactants
                    )

                    print(f"Is known sulfonamide reaction: {is_sulfonamide_rxn}")
                    print(f"Sulfonyl halide in reactants: {sulfonyl_halide_in_reactants}")
                    print(f"Sulfonamide in reactants: {sulfonamide_in_reactants}")
                    print(f"Amine in reactants: {amine_in_reactants}")

                    # Case 1: Classic sulfonamide formation (sulfonyl halide + amine)
                    if sulfonyl_halide_in_reactants and amine_in_reactants:
                        print(f"Confirmed late-stage sulfonamide formation at depth {depth}")
                        found_coupling = True

                    # Case 2: Coupling with a pre-existing sulfonamide group
                    elif sulfonamide_in_reactants:
                        print(
                            f"Confirmed late-stage coupling with sulfonamide fragment at depth {depth}"
                        )
                        found_coupling = True

                    # Case 3: Known sulfonamide reaction type
                    elif is_sulfonamide_rxn:
                        print(
                            f"Confirmed late-stage sulfonamide reaction (by reaction type) at depth {depth}"
                        )
                        found_coupling = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_coupling
