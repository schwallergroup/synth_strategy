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
    This function detects the use of ketone protection via ketal formation in the synthetic route.
    """
    ketone_protected = False

    def dfs_traverse(node):
        nonlocal ketone_protected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ketone protection (ketone + diol -> ketal)
            # Using the acetalization reaction check
            if checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ) or checker.check_reaction("Diol acetalization", rsmi):
                # Verify ketone in reactants
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants_smiles)

                # Check if product has a ketal structure
                has_ketal_product = checker.check_fg("Acetal/Ketal", product_smiles)

                if has_ketone and has_ketal_product:
                    print(f"Detected ketone protection via ketal formation: {rsmi}")
                    ketone_protected = True

            # Also check for ketal hydrolysis to ketone (deprotection)
            elif checker.check_reaction(
                "Ketal hydrolysis to ketone", rsmi
            ) or checker.check_reaction("Acetal hydrolysis to ketone", rsmi):
                # Verify ketal in reactants
                has_ketal = any(checker.check_fg("Acetal/Ketal", r) for r in reactants_smiles)

                # Check if product has a ketone structure
                has_ketone_product = checker.check_fg("Ketone", product_smiles)

                if has_ketal and has_ketone_product:
                    print(f"Detected ketone deprotection from ketal: {rsmi}")
                    ketone_protected = True

            # Alternative check if reaction type check fails
            elif not ketone_protected:
                # Check for ketone in reactants
                ketone_reactants = [r for r in reactants_smiles if checker.check_fg("Ketone", r)]

                # Check for ketal in product
                has_ketal_product = checker.check_fg("Acetal/Ketal", product_smiles)

                # Check for diol in reactants (two OH groups)
                diol_reactants = []
                for r in reactants_smiles:
                    # Check for alcohol functional groups
                    alcohol_count = 0
                    for fg in ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol"]:
                        if checker.check_fg(fg, r):
                            alcohol_indices = checker.get_fg_atom_indices(fg, r)
                            alcohol_count += len(alcohol_indices)

                    if alcohol_count >= 2:
                        diol_reactants.append(r)

                if ketone_reactants and diol_reactants and has_ketal_product:
                    print(f"Detected potential ketone protection: {rsmi}")
                    ketone_protected = True

                # Check for ketal in reactants and ketone in product (deprotection)
                ketal_reactants = [
                    r for r in reactants_smiles if checker.check_fg("Acetal/Ketal", r)
                ]
                has_ketone_product = checker.check_fg("Ketone", product_smiles)

                if ketal_reactants and has_ketone_product:
                    print(f"Detected potential ketone deprotection: {rsmi}")
                    ketone_protected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ketone_protected
