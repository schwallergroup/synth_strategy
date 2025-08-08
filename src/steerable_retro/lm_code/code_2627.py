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
    Detects if the synthesis route employs a protection-deprotection sequence,
    specifically looking for alcohol, amine, or carboxylic acid protection and deprotection.
    """
    # Track protection and deprotection steps with their depths and functional groups
    protection_steps = []
    deprotection_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection reactions using reaction checkers
            if (
                checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                or checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Protection of carboxylic acid", rsmi)
            ):
                # Identify which functional group is being protected
                fg_type = None
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    for r in reactants
                ):
                    fg_type = "alcohol"
                elif any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    fg_type = "amine"
                elif any(checker.check_fg("Carboxylic acid", r) for r in reactants):
                    fg_type = "carboxylic_acid"

                if fg_type:
                    protection_steps.append((depth, fg_type))
                    print(f"Found {fg_type} protection step at depth {depth}: {rsmi}")

            # Fallback to functional group checking for protection
            elif (
                any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    for r in reactants
                )
            ) and (
                checker.check_fg("Ester", product)
                or checker.check_fg("Silyl ether protective group", product)
                or checker.check_fg("TMS ether protective group", product)
            ):
                protection_steps.append((depth, "alcohol"))
                print(f"Found alcohol protection step (FG check) at depth {depth}: {rsmi}")

            elif (
                any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )
            ) and (checker.check_fg("Boc", product) or checker.check_fg("Acetal/Ketal", product)):
                protection_steps.append((depth, "amine"))
                print(f"Found amine protection step (FG check) at depth {depth}: {rsmi}")

            elif (any(checker.check_fg("Carboxylic acid", r) for r in reactants)) and (
                checker.check_fg("Ester", product)
            ):
                protection_steps.append((depth, "carboxylic_acid"))
                print(f"Found carboxylic acid protection step (FG check) at depth {depth}: {rsmi}")

            # Check for deprotection reactions using reaction checkers
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                or checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
                or checker.check_reaction("COOH ethyl deprotection", rsmi)
                or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
            ):

                # Identify which functional group is being deprotected
                fg_type = None
                if (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Phenol", product)
                ):
                    fg_type = "alcohol"
                elif checker.check_fg("Primary amine", product) or checker.check_fg(
                    "Secondary amine", product
                ):
                    fg_type = "amine"
                elif checker.check_fg("Carboxylic acid", product):
                    fg_type = "carboxylic_acid"

                if fg_type:
                    deprotection_steps.append((depth, fg_type))
                    print(f"Found {fg_type} deprotection step at depth {depth}: {rsmi}")

            # Fallback to functional group checking for deprotection
            elif (
                any(
                    checker.check_fg("Ester", r)
                    or checker.check_fg("Silyl ether protective group", r)
                    or checker.check_fg("TMS ether protective group", r)
                    for r in reactants
                )
            ) and (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
                or checker.check_fg("Phenol", product)
            ):
                deprotection_steps.append((depth, "alcohol"))
                print(f"Found alcohol deprotection step (FG check) at depth {depth}: {rsmi}")

            elif (any(checker.check_fg("Boc", r) for r in reactants)) and (
                checker.check_fg("Primary amine", product)
                or checker.check_fg("Secondary amine", product)
            ):
                deprotection_steps.append((depth, "amine"))
                print(f"Found amine deprotection step (FG check) at depth {depth}: {rsmi}")

            elif (any(checker.check_fg("Ester", r) for r in reactants)) and (
                checker.check_fg("Carboxylic acid", product)
            ):
                deprotection_steps.append((depth, "carboxylic_acid"))
                print(
                    f"Found carboxylic acid deprotection step (FG check) at depth {depth}: {rsmi}"
                )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both protection and deprotection steps
    if protection_steps and deprotection_steps:
        print(f"Protection steps at depths: {protection_steps}")
        print(f"Deprotection steps at depths: {deprotection_steps}")

        # Check for matching protection-deprotection pairs
        for prot_depth, prot_fg in protection_steps:
            for deprot_depth, deprot_fg in deprotection_steps:
                # In retrosynthesis, higher depth = earlier stage in forward synthesis
                # So for a proper protection-deprotection sequence, the protection depth
                # should be greater than the deprotection depth
                if prot_depth > deprot_depth and prot_fg == deprot_fg:
                    print(
                        f"Found {prot_fg} protection-deprotection sequence: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                    )
                    return True

    return False
