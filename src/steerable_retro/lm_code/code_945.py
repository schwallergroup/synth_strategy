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
    This function detects the use of silyl protection for alcohols.
    """
    found_silyl_ether = False
    protection_reaction = False

    def dfs_traverse(node):
        nonlocal found_silyl_ether, protection_reaction

        if node["type"] == "mol":
            # Check if the molecule contains a silyl ether protective group
            if checker.check_fg("TMS ether protective group", node["smiles"]) or checker.check_fg(
                "Silyl protective group", node["smiles"]
            ):
                found_silyl_ether = True
                print(f"Found silyl ether in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for silyl protection or deprotection reactions
                if (
                    checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rsmi
                    )
                    or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                ):
                    protection_reaction = True
                    print(f"Found silyl protection/deprotection reaction: {rsmi}")

                # If no explicit reaction type found, check for silyl protection pattern
                if not protection_reaction:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant is a silyl chloride
                        silyl_chloride_reactant = any(
                            checker.check_fg("Silyl protective group", r) and "Cl" in r
                            for r in reactants
                        )

                        # Check if any reactant has an alcohol and product has a silyl ether
                        alcohol_to_silyl = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            for r in reactants
                        ) and (
                            checker.check_fg("TMS ether protective group", product)
                            or checker.check_fg("Silyl protective group", product)
                        )

                        if silyl_chloride_reactant and alcohol_to_silyl:
                            protection_reaction = True
                            print(f"Found implicit silyl protection reaction: {rsmi}")

                        # Check for deprotection pattern (silyl ether to alcohol)
                        silyl_to_alcohol = any(
                            checker.check_fg("TMS ether protective group", r)
                            or checker.check_fg("Silyl protective group", r)
                            for r in reactants
                        ) and (
                            checker.check_fg("Primary alcohol", product)
                            or checker.check_fg("Secondary alcohol", product)
                            or checker.check_fg("Tertiary alcohol", product)
                        )

                        if silyl_to_alcohol:
                            protection_reaction = True
                            print(f"Found implicit silyl deprotection reaction: {rsmi}")

                    except Exception as e:
                        print(f"Error analyzing reaction components: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = found_silyl_ether and protection_reaction
    print(
        f"Found silyl ether: {found_silyl_ether}, Found protection reaction: {protection_reaction}"
    )
    return result
