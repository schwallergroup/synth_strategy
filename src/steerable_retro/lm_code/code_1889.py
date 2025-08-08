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
    This function detects if the synthetic route employs aromatic cyanation,
    specifically looking for installation of a cyano group on an aryl halide.
    """
    cyanation_found = False

    def dfs_traverse(node):
        nonlocal cyanation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide in reactants (Br, Cl, I, F)
            aryl_halide_present = any(checker.check_fg("Aromatic halide", r) for r in reactants)

            # Check for nitrile group in product
            nitrile_in_product = checker.check_fg("Nitrile", product)

            # Check for cyanide source in reactants (optional, as it might be implicit)
            cyanide_source = any(checker.check_fg("Nitrile", r) for r in reactants)

            # Verify this is a cyanation reaction
            if aryl_halide_present and nitrile_in_product:
                # If there's a cyanide source in reactants and it's not the main reactant with the aryl halide
                if cyanide_source:
                    # Make sure at least one reactant has a nitrile but no aromatic halide
                    # (this would be the cyanide source)
                    has_cyanide_source = False
                    for r in reactants:
                        if checker.check_fg("Nitrile", r) and not checker.check_fg(
                            "Aromatic halide", r
                        ):
                            has_cyanide_source = True
                            break

                    if has_cyanide_source:
                        print(
                            f"Found aromatic cyanation reaction with explicit cyanide source: {rsmi}"
                        )
                        cyanation_found = True
                else:
                    # Even without an explicit cyanide source, if we have an aryl halide in reactants
                    # and a nitrile in the product, it's likely a cyanation reaction
                    print(
                        f"Found potential aromatic cyanation reaction without explicit cyanide source: {rsmi}"
                    )
                    cyanation_found = True

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return cyanation_found
