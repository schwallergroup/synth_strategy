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
    Detects if the synthesis includes an amide disconnection step.
    """
    has_amide_disconnection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_disconnection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if this is an amide hydrolysis reaction
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                if (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):
                    print(f"Found amide disconnection via hydrolysis at depth {depth}")
                    has_amide_disconnection = True
                    return

            # Check for amide formation reactions
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
            ):
                print(f"Found amide formation reaction at depth {depth}")
                has_amide_disconnection = True
                return

            # Manual check for amide disconnection if specific reaction types not detected
            # In retrosynthesis: product (amide) -> reactants (amine + carboxylic acid or derivative)
            if (
                checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Secondary amide", product_smiles)
                or checker.check_fg("Tertiary amide", product_smiles)
            ):

                # Check if reactants contain amine and carboxylic acid or derivative
                has_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants_smiles.split(".")
                )

                has_carboxylic = any(
                    checker.check_fg("Carboxylic acid", r)
                    or checker.check_fg("Acyl halide", r)
                    or checker.check_fg("Ester", r)
                    or checker.check_fg("Anhydride", r)
                    for r in reactants_smiles.split(".")
                )

                if has_amine and has_carboxylic:
                    print(f"Found amide disconnection via manual check at depth {depth}")
                    has_amide_disconnection = True
                    return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_amide_disconnection
