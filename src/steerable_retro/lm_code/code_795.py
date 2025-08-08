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
    This function detects a strategy involving sulfonyl group coupling.
    """
    found_sulfonyl_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonyl_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Direct check for known sulfonyl coupling reactions
            if (
                checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )
                or checker.check_reaction("Formation of Sulfonic Esters", rsmi)
                or checker.check_reaction(
                    "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
                )
            ):
                print(f"Found sulfonyl coupling reaction at depth {depth}: {rsmi}")
                found_sulfonyl_coupling = True
                return

            # If no direct reaction match, check for structural changes
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonyl groups in reactants
            sulfonyl_in_reactants = False
            for r in reactants:
                if (
                    checker.check_fg("Sulfonamide", r)
                    or checker.check_fg("Sulfone", r)
                    or checker.check_fg("Sulfonate", r)
                    or checker.check_fg("Sulfonic acid", r)
                    or checker.check_fg("Sulfonyl halide", r)
                ):
                    sulfonyl_in_reactants = True
                    break

            # If sulfonyl group in reactants, check for new C-S bond formation in product
            if sulfonyl_in_reactants:
                # Check if product has a sulfonyl-containing group
                if (
                    checker.check_fg("Sulfonamide", product)
                    or checker.check_fg("Sulfone", product)
                    or checker.check_fg("Sulfonate", product)
                ):

                    # Check if this is likely a coupling reaction
                    # Look for reactants that could be coupling partners (e.g., alcohols, amines)
                    potential_coupling_partner = False
                    for r in reactants:
                        if (
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Phenol", r)
                        ):
                            potential_coupling_partner = True
                            break

                    if potential_coupling_partner:
                        print(f"Found sulfonyl coupling at depth {depth}: {rsmi}")
                        found_sulfonyl_coupling = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_sulfonyl_coupling
