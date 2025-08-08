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
    Detects a strategy involving late-stage esterification of a carboxylic acid
    after amide formation with a cyano-containing fragment.
    """
    # Track if we found the key features
    found_esterification = False
    found_amide_formation = False
    found_cyano_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal found_esterification, found_amide_formation, found_cyano_fragment

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            if not product_smiles or not reactants_smiles:
                return

            # Check for esterification at depth 0-3 (late stage)
            if depth <= 3:
                # Look for carboxylic acid in reactants and ester in product
                has_acid = any(
                    checker.check_fg("Carboxylic acid", r_smiles)
                    for r_smiles in reactants_smiles
                    if r_smiles
                )
                has_alcohol = any(
                    checker.check_fg(alcohol_type, r_smiles)
                    for r_smiles in reactants_smiles
                    if r_smiles
                    for alcohol_type in [
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Aromatic alcohol",
                    ]
                )
                has_ester = checker.check_fg("Ester", product_smiles) if product_smiles else False
                has_acyl_halide = any(
                    checker.check_fg("Acyl halide", r_smiles)
                    for r_smiles in reactants_smiles
                    if r_smiles
                )

                # Check for various esterification reactions
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Transesterification",
                    "Oxidative esterification of primary alcohols",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Acetic anhydride and alcohol to ester",
                ]

                is_esterification = any(
                    checker.check_reaction(rxn, rsmi) for rxn in esterification_reactions
                )

                # Check for conversion of acid to ester without specific reaction type
                acid_to_ester = has_acid and has_ester

                # Check for acyl halide to ester conversion
                acyl_to_ester = has_acyl_halide and has_ester and has_alcohol

                # Also check for direct ester formation
                direct_ester_formation = (
                    not has_acid
                    and has_ester
                    and not any(
                        checker.check_fg("Ester", r_smiles)
                        for r_smiles in reactants_smiles
                        if r_smiles
                    )
                )

                if is_esterification or acid_to_ester or acyl_to_ester or direct_ester_formation:
                    found_esterification = True
                    print(f"Found late-stage esterification at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")

            # Check for amide formation with cyano fragment
            if depth >= 1:
                # Check for amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "{Schotten-Baumann_amide}",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in amide_formation_reactions
                )

                # Check if product has amide bond
                has_amide = any(
                    checker.check_fg(amide_type, product_smiles)
                    for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]
                )

                # Check if any reactant has cyano group or if product has both amide and cyano
                has_cyano_reactant = any(
                    checker.check_fg("Nitrile", r_smiles)
                    for r_smiles in reactants_smiles
                    if r_smiles
                )
                has_cyano_product = (
                    checker.check_fg("Nitrile", product_smiles) if product_smiles else False
                )

                if is_amide_formation or has_amide:
                    found_amide_formation = True
                    print(f"Found amide formation at depth {depth}")

                if has_cyano_reactant or has_cyano_product:
                    found_cyano_fragment = True
                    print(f"Found cyano-containing fragment at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if all key features were found
    result = found_esterification and found_amide_formation and found_cyano_fragment
    print(f"Strategy detection result: {result}")
    print(
        f"Esterification: {found_esterification}, Amide formation: {found_amide_formation}, Cyano fragment: {found_cyano_fragment}"
    )

    return result
