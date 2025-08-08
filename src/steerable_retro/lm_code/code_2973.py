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
    This function detects the use of TBDMS as a protecting group for alcohols,
    including both protection and deprotection steps.
    """
    found_tbdms_group = False
    found_protection = False
    found_deprotection = False

    def dfs_traverse(node):
        nonlocal found_tbdms_group, found_protection, found_deprotection

        if node["type"] == "mol":
            # Check if molecule contains silyl-protected alcohol (TBDMS is a type of silyl group)
            if checker.check_fg("Silyl protective group", node["smiles"]):
                # Additional check to confirm it's specifically TBDMS
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and (
                    "[Si](C)(C)C(C)(C)C" in node["smiles"]
                    or "C[Si](C)C(C)(C)C" in node["smiles"]
                    or "[Si]([C])([C])C([C])([C])[C]" in node["smiles"]
                ):
                    found_tbdms_group = True
                    print(f"Confirmed TBDMS structure in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for TBDMS protection
                if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    # Verify an alcohol is being protected
                    alcohol_present = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                        ):
                            alcohol_present = True
                            break

                    # Verify TBDMS reagent is present
                    tbdms_reagent_present = False
                    for reactant in reactants:
                        # Look for silicon-containing reagent with t-butyl group
                        if (
                            "[Si]" in reactant
                            and "Cl" in reactant
                            and ("C(C)(C)C" in reactant or "CC(C)(C)" in reactant)
                        ):
                            tbdms_reagent_present = True
                            break

                    # Verify product contains TBDMS group
                    tbdms_product = checker.check_fg("Silyl protective group", product) and (
                        "[Si](C)(C)C(C)(C)C" in product
                        or "C[Si](C)C(C)(C)C" in product
                        or "[Si]([C])([C])C([C])([C])[C]" in product
                    )

                    if alcohol_present and tbdms_reagent_present and tbdms_product:
                        found_protection = True
                        print(f"Found TBDMS protection step: {rsmi}")

                # Check for TBDMS deprotection
                if (
                    checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rsmi
                    )
                    or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                ):

                    # Verify a TBDMS-protected alcohol is being deprotected
                    tbdms_reactant = False
                    for reactant in reactants:
                        if checker.check_fg("Silyl protective group", reactant) and (
                            "[Si](C)(C)C(C)(C)C" in reactant
                            or "C[Si](C)C(C)(C)C" in reactant
                            or "[Si]([C])([C])C([C])([C])[C]" in reactant
                        ):
                            tbdms_reactant = True
                            break

                    # Verify product is an alcohol
                    alcohol_product = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )

                    if tbdms_reactant and alcohol_product:
                        found_deprotection = True
                        print(f"Found TBDMS deprotection step: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if TBDMS group is found along with either protection or deprotection step
    result = found_tbdms_group and (found_protection or found_deprotection)
    print(
        f"TBDMS group found: {found_tbdms_group}, Protection found: {found_protection}, Deprotection found: {found_deprotection}"
    )

    # Additional check: if we found TBDMS-Cl as a molecule, it's likely part of a protection strategy
    if not result and found_tbdms_group:
        # Check if any molecule is TBDMS-Cl
        def check_for_tbdms_cl(node):
            if node["type"] == "mol" and "Cl" in node["smiles"] and "[Si]" in node["smiles"]:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and ("C(C)(C)C" in node["smiles"] or "CC(C)(C)" in node["smiles"]):
                    print(f"Found TBDMS-Cl reagent: {node['smiles']}")
                    return True

            for child in node.get("children", []):
                if check_for_tbdms_cl(child):
                    return True
            return False

        if check_for_tbdms_cl(route):
            print("TBDMS protection strategy detected based on presence of TBDMS-Cl reagent")
            result = True

    return result
