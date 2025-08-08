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
    This function detects a strategy involving amide coupling between a heterocyclic
    carboxylic acid and an amine with a stereocenter.
    """
    # Define heterocyclic rings to check
    heterocyclic_rings = [
        "furan",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
    ]

    # Track if we've found the required components in a single reaction
    found_valid_reaction = False

    def dfs_traverse(node):
        nonlocal found_valid_reaction

        if node["type"] == "reaction" and not found_valid_reaction:
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling reaction
            is_amide_coupling = (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Schotten-Baumann_amide", rsmi)
            )

            if not is_amide_coupling:
                # Alternative check for amide formation
                has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                has_amine = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants_smiles
                )
                has_amide = checker.check_fg("Primary amide", product_smiles) or checker.check_fg(
                    "Secondary amide", product_smiles
                )

                is_amide_coupling = has_acid and has_amine and has_amide

            if is_amide_coupling:
                print(f"Found amide coupling reaction: {rsmi}")

                # Check for heterocyclic acid
                heterocyclic_acid_found = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", reactant):
                        # Check if it contains a heterocyclic ring
                        for ring in heterocyclic_rings:
                            if checker.check_ring(ring, reactant):
                                print(f"Found heterocyclic acid with {ring} ring: {reactant}")
                                heterocyclic_acid_found = True
                                break

                # Check for amine with stereocenter
                stereocenter_amine_found = False
                for reactant in reactants_smiles:
                    if ("@H" in reactant or "@@H" in reactant) and (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                    ):
                        print(f"Found amine with stereocenter: {reactant}")
                        stereocenter_amine_found = True

                # If both conditions are met in the same reaction
                if heterocyclic_acid_found and stereocenter_amine_found:
                    print("Found valid heterocyclic amide coupling with stereocenter")
                    found_valid_reaction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_valid_reaction
