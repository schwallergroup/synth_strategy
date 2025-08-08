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
    This function detects an ester hydrolysis followed by amide formation strategy,
    where an ester is first hydrolyzed to an acid and then converted to an amide.
    """
    found_ester_hydrolysis = False
    found_amide_formation = False
    hydrolysis_depth = -1
    amide_depth = -1
    acid_from_hydrolysis = None

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis, found_amide_formation, hydrolysis_depth, amide_depth, acid_from_hydrolysis

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for ester hydrolysis - expanded to include more reaction types
                hydrolysis_reactions = [
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "COOH ethyl deprotection",
                ]

                is_hydrolysis = any(
                    checker.check_reaction(rxn, rsmi) for rxn in hydrolysis_reactions
                )

                # If reaction check fails, try to detect by functional group transformation
                if not is_hydrolysis:
                    reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants if r)
                    product_has_acid = checker.check_fg("Carboxylic acid", product)
                    if (
                        reactant_has_ester
                        and product_has_acid
                        and not any(checker.check_fg("Carboxylic acid", r) for r in reactants if r)
                    ):
                        print(f"Detected ester hydrolysis by functional group transformation")
                        is_hydrolysis = True

                if is_hydrolysis:
                    print(f"Found ester hydrolysis at depth {depth}")
                    # Verify reactant is an ester and product is an acid
                    reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants if r)
                    product_has_acid = checker.check_fg("Carboxylic acid", product)

                    if reactant_has_ester and product_has_acid:
                        print(f"Confirmed ester hydrolysis: ester → acid")
                        found_ester_hydrolysis = True
                        hydrolysis_depth = depth
                        # Store the acid product for later verification
                        acid_from_hydrolysis = product
                        print(f"Stored acid from hydrolysis: {acid_from_hydrolysis}")

                # Check for amide formation reactions - expanded list
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Carboxylic acid to amide conversion",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in amide_formation_reactions
                )

                # If reaction check fails, try to detect by functional group transformation
                if not is_amide_formation:
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants if r)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                        if r
                    )
                    has_amide_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if has_acid and has_amine and has_amide_product:
                        print(f"Detected amide formation by functional group transformation")
                        is_amide_formation = True

                if is_amide_formation:
                    print(f"Found potential amide formation at depth {depth}")
                    # Check reactants and products
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants if r)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                        if r
                    )
                    has_amide_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    print(
                        f"Amide formation check - Has acid: {has_acid}, Has amine: {has_amine}, Has amide product: {has_amide_product}"
                    )

                    if (
                        (has_acid or any(checker.check_fg("Ester", r) for r in reactants if r))
                        and has_amine
                        and has_amide_product
                    ):
                        print(f"Confirmed amide formation: acid/ester + amine → amide")
                        found_amide_formation = True
                        amide_depth = depth

                        # If we found hydrolysis earlier, we've confirmed the pattern
                        if found_ester_hydrolysis:
                            print(f"Found both hydrolysis and amide formation")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # In retrosynthetic traversal, later steps have lower depth values
    # So amide formation should have a lower depth than hydrolysis
    if found_ester_hydrolysis and found_amide_formation:
        print(f"Hydrolysis depth: {hydrolysis_depth}, Amide formation depth: {amide_depth}")
        return amide_depth < hydrolysis_depth

    print(
        f"Did not find the pattern. Ester hydrolysis: {found_ester_hydrolysis}, Amide formation: {found_amide_formation}"
    )
    return False
