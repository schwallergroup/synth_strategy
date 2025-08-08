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
    This function detects if the synthesis route uses THP protection and deprotection
    for alcohol protection.
    """
    has_thp_protection = False
    has_thp_deprotection = False

    def dfs_traverse(node):
        nonlocal has_thp_protection, has_thp_deprotection

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for THP protection (alcohol + THP reagent -> THP protected alcohol)
                if not has_thp_protection:
                    # Check for alcohol in reactants
                    has_alcohol = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            has_alcohol = True
                            break

                    # Check for THP ring or dihydropyran in reactants
                    has_thp_reagent = False
                    for reactant in reactants:
                        # Check for tetrahydropyran ring or dihydropyran (3,4-dihydro-2H-pyran)
                        if (
                            checker.check_ring("tetrahydropyran", reactant)
                            or "O1CCCC=C1" in reactant
                        ):
                            has_thp_reagent = True
                            break

                    # Check for THP ether in product
                    has_thp_product = checker.check_ring("tetrahydropyran", product)

                    # Check if this is an alcohol protection reaction
                    if has_alcohol and has_thp_reagent and has_thp_product:
                        print(f"THP protection found: {rsmi}")
                        has_thp_protection = True

                # Check for THP deprotection (THP protected alcohol -> alcohol)
                if not has_thp_deprotection:
                    # Check for THP ring in reactants
                    has_thp_reactant = False
                    for reactant in reactants:
                        if checker.check_ring("tetrahydropyran", reactant):
                            has_thp_reactant = True
                            break

                    # Check for alcohol in product
                    has_alcohol_product = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    )

                    # Check if this is an alcohol deprotection reaction
                    if has_thp_reactant and has_alcohol_product:
                        print(f"THP deprotection found: {rsmi}")
                        has_thp_deprotection = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Protection found: {has_thp_protection}, Deprotection found: {has_thp_deprotection}")

    # In retrosynthesis, we might only see one direction of the protection/deprotection
    # Either one indicates THP was used in the synthesis
    return has_thp_protection or has_thp_deprotection
