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
    Detects if the synthesis involves multiple amide bond formations.
    """
    amide_formation_count = 0
    amide_formation_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Carboxylic acid with primary amine to amide",
        "Acyl chloride with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Ester with ammonia to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
        "Acylation of primary amines",
        "Acylation of secondary amines",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an amide formation reaction by reaction type
                is_amide_formation = False
                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type}")
                        is_amide_formation = True
                        break

                # If not identified by reaction type, check by functional group changes
                if not is_amide_formation:
                    # Count amides in product
                    product_amide_count = 0
                    if checker.check_fg("Primary amide", product):
                        product_amide_count += 1
                    if checker.check_fg("Secondary amide", product):
                        product_amide_count += 1
                    if checker.check_fg("Tertiary amide", product):
                        product_amide_count += 1

                    # Count amides in reactants
                    reactant_amide_count = 0
                    for reactant in reactants:
                        if checker.check_fg("Primary amide", reactant):
                            reactant_amide_count += 1
                        if checker.check_fg("Secondary amide", reactant):
                            reactant_amide_count += 1
                        if checker.check_fg("Tertiary amide", reactant):
                            reactant_amide_count += 1

                    # Check if an amide was formed
                    if product_amide_count > reactant_amide_count:
                        print(
                            f"Detected amide formation by FG count: product={product_amide_count}, reactants={reactant_amide_count}"
                        )

                        # Verify this is actually an amide formation by checking for reactants
                        has_amine = False
                        has_carboxylic_acid_or_derivative = False

                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Tertiary amine", reactant)
                            ):
                                has_amine = True

                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Anhydride", reactant)
                            ):
                                has_carboxylic_acid_or_derivative = True

                        if has_amine and has_carboxylic_acid_or_derivative:
                            is_amide_formation = True
                            print(f"Confirmed amide formation by reactant analysis")

                if is_amide_formation:
                    amide_formation_count += 1
                    print(f"Total amide formations so far: {amide_formation_count}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final amide formation count: {amide_formation_count}")
    return amide_formation_count >= 2
