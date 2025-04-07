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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthesis involves amide formation from an acid chloride
    and an amine.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if the reaction is an acylation of nitrogen nucleophiles by acyl halides
            if checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                rsmi,
            ):
                print(
                    "Found amide formation from acid chloride and amine using reaction type check"
                )
                amide_formation_found = True
            else:
                # Fallback to manual checking if reaction type check fails
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride in reactants
                acid_chloride_present = any(
                    checker.check_fg("Acyl halide", reactant) for reactant in reactants
                )

                # Check for amine in reactants
                amine_present = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    for reactant in reactants
                )

                # Check for amide in product
                amide_present = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                # Check for Schotten-Baumann reaction (another way to form amides from acid chlorides)
                schotten_baumann = (
                    checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                )

                if (acid_chloride_present and amine_present and amide_present) or schotten_baumann:
                    print(
                        "Found amide formation from acid chloride and amine using functional group checks"
                    )
                    amide_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
