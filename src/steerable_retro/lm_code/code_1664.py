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
    This function detects esterification of carboxylic acid or related transformations
    including ester hydrolysis, saponification, and transesterification.
    """
    has_esterification = False

    def dfs_traverse(node):
        nonlocal has_esterification

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction: {rsmi}")

                # Check for direct esterification
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    print(f"Detected esterification reaction: {rsmi}")
                    has_esterification = True
                    return

                # Check for ester hydrolysis (reverse of esterification)
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    print(f"Detected ester hydrolysis reaction: {rsmi}")
                    has_esterification = True
                    return

                # Check for ester saponification
                if checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                    print(f"Detected ester saponification reaction: {rsmi}")
                    has_esterification = True
                    return

                # Check for transesterification
                if checker.check_reaction("Transesterification", rsmi):
                    print(f"Detected transesterification reaction: {rsmi}")
                    has_esterification = True
                    return

                # Check for O-alkylation with diazo compounds
                if checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                ):
                    print(f"Detected diazo esterification: {rsmi}")
                    has_esterification = True
                    return

                print("Not identified by reaction checker, trying functional group analysis")

                # Check for carboxylic acid in reactants and ester in product (esterification)
                has_carboxylic_acid_reactants = any(
                    checker.check_fg("Carboxylic acid", reactant) for reactant in reactants
                )
                has_ester_product = checker.check_fg("Ester", product)

                # Check for ester in reactants and carboxylic acid in product (hydrolysis)
                has_ester_reactants = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants
                )
                has_carboxylic_acid_product = checker.check_fg("Carboxylic acid", product)

                # Check for alcohol in reactants (typical for esterification)
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", reactant)
                    or checker.check_fg("Secondary alcohol", reactant)
                    or checker.check_fg("Tertiary alcohol", reactant)
                    or checker.check_fg("Aromatic alcohol", reactant)
                    or checker.check_fg("Enol", reactant)
                    for reactant in reactants
                )

                # Check for diazo compounds (alternative esterification route)
                has_diazo = any(checker.check_fg("Diazo", reactant) for reactant in reactants)

                print(f"Has carboxylic acid in reactants: {has_carboxylic_acid_reactants}")
                print(f"Has ester in product: {has_ester_product}")
                print(f"Has ester in reactants: {has_ester_reactants}")
                print(f"Has carboxylic acid in product: {has_carboxylic_acid_product}")
                print(f"Has alcohol in reactants: {has_alcohol}")
                print(f"Has diazo compound in reactants: {has_diazo}")

                # Esterification: carboxylic acid + alcohol → ester
                if has_carboxylic_acid_reactants and has_ester_product:
                    print(
                        f"Detected functional group changes consistent with esterification: {rsmi}"
                    )
                    has_esterification = True
                    return

                # Hydrolysis: ester → carboxylic acid
                if has_ester_reactants and has_carboxylic_acid_product:
                    print(
                        f"Detected functional group changes consistent with ester hydrolysis: {rsmi}"
                    )
                    has_esterification = True
                    return

                # Transesterification: ester + alcohol → different ester
                if has_ester_reactants and has_ester_product and has_alcohol:
                    print(
                        f"Detected functional group changes consistent with transesterification: {rsmi}"
                    )
                    has_esterification = True
                    return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_esterification
