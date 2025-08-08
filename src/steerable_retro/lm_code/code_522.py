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
    This function detects the formation of an isoxazoline ring via 1,3-dipolar cycloaddition
    between an oxime and an alkene.
    """
    found_cycloaddition = False

    def dfs_traverse(node):
        nonlocal found_cycloaddition

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check if reactants contain oxime and alkene
                reactants_have_oxime = any(checker.check_fg("Oxime", r) for r in reactants_smiles)
                reactants_have_alkene = any(
                    checker.check_fg("Alkene", r)
                    or checker.check_fg("Vinyl", r)
                    or checker.check_fg("Ethylene", r)
                    for r in reactants_smiles
                )

                print(f"Reactants have oxime: {reactants_have_oxime}")
                print(f"Reactants have alkene: {reactants_have_alkene}")

                # Check for isoxazoline structure in product
                product_mol = Chem.MolFromSmiles(product_smiles)
                product_smiles_canonical = Chem.MolToSmiles(product_mol)

                # Check if this is a 1,3-dipolar cycloaddition reaction
                is_cycloaddition = (
                    checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                    or checker.check_reaction(
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction("[3+2]-cycloaddition of hydrazone and alkene", rsmi)
                )

                print(f"Is cycloaddition reaction: {is_cycloaddition}")

                # Check if oxime is consumed (not present in product)
                product_has_oxime = checker.check_fg("Oxime", product_smiles)
                print(f"Product has oxime: {product_has_oxime}")

                # Manual check for isoxazoline structure (5-membered ring with N-O bond)
                # Isoxazoline pattern: C1=NOC(C)C1 or similar patterns
                isoxazoline_pattern = Chem.MolFromSmarts("C1=NOC(C)C1")
                isoxazoline_pattern2 = Chem.MolFromSmarts("C1=NOC(*)C1")
                isoxazoline_pattern3 = Chem.MolFromSmarts("C1=NOC([*])([*])C1")

                has_isoxazoline_structure = (
                    product_mol.HasSubstructMatch(isoxazoline_pattern)
                    or product_mol.HasSubstructMatch(isoxazoline_pattern2)
                    or product_mol.HasSubstructMatch(isoxazoline_pattern3)
                )

                # Check for N-O bond in a ring
                has_no_in_ring = (
                    "N1OC" in product_smiles_canonical or "NOC1" in product_smiles_canonical
                )

                print(f"Has isoxazoline structure: {has_isoxazoline_structure}")
                print(f"Has N-O bond in ring: {has_no_in_ring}")

                # Criteria for isoxazoline formation:
                # 1. Reactants have oxime and alkene
                # 2. Product has isoxazoline structure
                # 3. Either it's a known cycloaddition reaction OR oxime is consumed
                if (
                    reactants_have_oxime
                    and reactants_have_alkene
                    and (has_isoxazoline_structure or has_no_in_ring)
                    and (is_cycloaddition or (reactants_have_oxime and not product_has_oxime))
                ):
                    print(f"Confirmed isoxazoline formation via 1,3-dipolar cycloaddition")
                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")
                    found_cycloaddition = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_cycloaddition
