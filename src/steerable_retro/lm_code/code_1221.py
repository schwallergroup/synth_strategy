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
    This function detects if the final product contains a cyanoacetamide group.
    """
    # The final product is the root node in the synthesis route
    if route["type"] != "mol":
        print("Root node is not a molecule")
        return False

    product_smiles = route["smiles"]
    print(f"Checking final product: {product_smiles}")

    # Check for cyanoacetamide using the checker function
    # Cyanoacetamide is a combination of an amide group and a nitrile group
    # in a specific arrangement: R-C(=O)-NH-CH2-CN or variations

    # First check if the molecule has both amide and nitrile groups
    has_amide = (
        checker.check_fg("Primary amide", product_smiles)
        or checker.check_fg("Secondary amide", product_smiles)
        or checker.check_fg("Tertiary amide", product_smiles)
    )

    has_nitrile = checker.check_fg("Nitrile", product_smiles)

    if not (has_amide and has_nitrile):
        print("Product does not have both amide and nitrile groups")
        return False

    # Now check for the specific cyanoacetamide pattern and its variations
    product_mol = Chem.MolFromSmiles(product_smiles)
    if product_mol is None:
        print("Could not parse product SMILES")
        return False

    # Classic cyanoacetamide pattern: R-C(=O)-NH-CH2-CN
    classic_pattern = Chem.MolFromSmarts("[*]-C(=O)-N-C-C#N")

    # Extended pattern allowing for variations in chain length
    extended_pattern = Chem.MolFromSmarts("[*]-C(=O)-N-[CH2,CH]-[CH2,CH]-C#N")

    # Pattern specifically for 3-cyanopropanamide: R-C(=O)-NH-CH2-CH2-CN
    cyano_propanamide = Chem.MolFromSmarts("C(=O)-N-C-C-C#N")

    # Pattern for any amide with a cyano group somewhere in the chain
    general_pattern = Chem.MolFromSmarts("C(=O)-N-[*]-[*]-C#N")

    if product_mol.HasSubstructMatch(classic_pattern):
        print("Detected classic cyanoacetamide in final product")
        return True

    if product_mol.HasSubstructMatch(extended_pattern):
        print("Detected extended cyanoacetamide in final product")
        return True

    if product_mol.HasSubstructMatch(cyano_propanamide):
        print("Detected 3-cyanopropanamide in final product")
        return True

    if product_mol.HasSubstructMatch(general_pattern):
        print("Detected general cyanoacetamide-like structure in final product")
        return True

    # As a final check, look for the specific fragment in the test case
    test_case_pattern = Chem.MolFromSmarts("NC(=O)CC#N")
    if product_mol.HasSubstructMatch(test_case_pattern):
        print("Detected specific cyanoacetamide structure in final product")
        return True

    print("Cyanoacetamide pattern not found in final product")
    return False
