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
    This function specifically detects the use of Boc protection strategy in the synthesis.
    It looks for Boc protection and deprotection events.
    """
    boc_protection_count = 0
    boc_deprotection_count = 0
    molecules_with_boc = set()  # Track molecules containing Boc groups

    # List of Boc protection and deprotection reaction types
    boc_protection_reactions = [
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
    ]

    boc_deprotection_reactions = [
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Tert-butyl deprotection of amine",
    ]

    # SMARTS pattern for Boc group (tert-butyloxycarbonyl)
    # We'll use functional group checking instead

    def check_molecule_for_boc(smiles):
        """Check if a molecule contains a Boc group"""
        # Check for carbamic ester (Boc is a type of carbamic ester)
        has_carbamic_ester = checker.check_fg("Carbamic ester", smiles)

        # Additional check for tert-butyl group connected to carbamic ester
        has_tertiary_alcohol = checker.check_fg("Tertiary alcohol", smiles)

        return has_carbamic_ester or has_tertiary_alcohol

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_count, boc_deprotection_count

        # Check molecule nodes for Boc groups
        if node["type"] == "mol" and node["smiles"]:
            if check_molecule_for_boc(node["smiles"]):
                molecules_with_boc.add(node["smiles"])
                print(f"Molecule with potential Boc group found: {node['smiles']}")

        # Check reaction nodes for Boc protection/deprotection
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection reactions
                for reaction_type in boc_protection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # Verify: product should have Boc group, at least one reactant shouldn't
                        if check_molecule_for_boc(product) and not all(
                            check_molecule_for_boc(r) for r in reactants
                        ):
                            boc_protection_count += 1
                            print(f"Boc protection detected ({reaction_type}): {rsmi}")
                            break

                # Check for Boc deprotection reactions
                for reaction_type in boc_deprotection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # Verify: at least one reactant should have Boc group, product shouldn't
                        if any(
                            check_molecule_for_boc(r) for r in reactants
                        ) and not check_molecule_for_boc(product):
                            boc_deprotection_count += 1
                            print(f"Boc deprotection detected ({reaction_type}): {rsmi}")
                            break
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have at least one Boc protection OR deprotection OR molecules with Boc
    has_boc_strategy = (
        boc_protection_count >= 1 or boc_deprotection_count >= 1 or len(molecules_with_boc) > 0
    )
    print(
        f"Boc protection count: {boc_protection_count}, Boc deprotection count: {boc_deprotection_count}, Molecules with Boc: {len(molecules_with_boc)}"
    )
    return has_boc_strategy
