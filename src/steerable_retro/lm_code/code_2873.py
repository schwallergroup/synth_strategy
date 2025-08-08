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
    This function detects amide bond formation via acyl chloride + amine reaction.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an acylation reaction using the checker function
                if checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                ):
                    print(f"Found acylation reaction: {rsmi}")
                    amide_formation_found = True
                    return

                # If the specific reaction check fails, try checking for the functional groups
                acyl_chloride_present = False
                amine_present = False
                amide_present = False

                # Check for amide in product
                if (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                ):
                    amide_present = True
                    print(f"Found amide in product: {product}")

                # Check for acyl chloride and amine in reactants
                for reactant in reactants:
                    if checker.check_fg("Acyl halide", reactant):
                        acyl_chloride_present = True
                        print(f"Found acyl halide in reactant: {reactant}")

                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        amine_present = True
                        print(f"Found amine in reactant: {reactant}")

                # Check if we have all the necessary components
                if amide_present and acyl_chloride_present and amine_present:
                    print(
                        "Amide formation from acyl chloride detected through functional group analysis"
                    )
                    amide_formation_found = True

                # Also check for specific reaction patterns
                if (
                    checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                ):
                    print(f"Found Schotten-Baumann or related reaction: {rsmi}")
                    amide_formation_found = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            if not amide_formation_found:  # Only continue if we haven't found the reaction yet
                dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
