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
    This function detects cyanation reactions (introduction of C≡N group).
    """
    cyanation_found = False

    def dfs_traverse(node):
        nonlocal cyanation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has nitrile group
            has_nitrile_in_product = checker.check_fg("Nitrile", product_smiles)

            if has_nitrile_in_product:
                # Count reactants with nitrile
                reactants_with_nitrile = sum(
                    1 for r_smiles in reactants_smiles if checker.check_fg("Nitrile", r_smiles)
                )

                # Check for common cyanation reagents
                cyanation_reagent_present = False
                for r_smiles in reactants_smiles:
                    # Common cyanation reagents
                    if any(
                        reagent in r_smiles.upper()
                        for reagent in ["CN", "KCN", "NACN", "CUCN", "TMSCN", "ZN(CN)2"]
                    ):
                        cyanation_reagent_present = True
                        break

                # Check if this is a cyanation reaction
                if cyanation_reagent_present and reactants_with_nitrile < len(reactants_smiles):
                    print(f"Cyanation detected with cyanation reagent: {rsmi}")
                    cyanation_found = True
                    return

                # Check if any reactant has a halide that could be replaced by CN
                for r_smiles in reactants_smiles:
                    if (
                        checker.check_fg("Primary halide", r_smiles)
                        or checker.check_fg("Secondary halide", r_smiles)
                        or checker.check_fg("Tertiary halide", r_smiles)
                        or checker.check_fg("Aromatic halide", r_smiles)
                    ):

                        # If reactant has halide and product has nitrile where reactant didn't
                        if not checker.check_fg("Nitrile", r_smiles):
                            print(f"Cyanation detected with halide replacement: {rsmi}")
                            cyanation_found = True
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return cyanation_found
