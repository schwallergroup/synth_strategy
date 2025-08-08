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
    Detects synthesis strategies that preserve key functional groups (cyano, methoxy)
    throughout the synthesis without modification.
    """
    # Track functional groups through synthesis
    cyano_preserved = True
    methoxy_preserved = True
    found_cyano = False
    found_methoxy = False

    def dfs_traverse(node, depth=0):
        nonlocal cyano_preserved, methoxy_preserved, found_cyano, found_methoxy

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants_smiles = reactants_part.split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for functional groups in reactants and product
            reactants_have_cyano = any(checker.check_fg("Nitrile", r) for r in reactants_smiles)
            product_has_cyano = checker.check_fg("Nitrile", product_smiles)

            # For methoxy, we need to check for ethers and specifically look for OCH3 groups
            reactants_have_methoxy = any(is_methoxy_present(r) for r in reactants_smiles)
            product_has_methoxy = is_methoxy_present(product_smiles)

            # Update preservation flags - since we're traversing retrosynthetically,
            # we need to check if groups in the product are preserved in the reactants
            if product_has_cyano and not reactants_have_cyano:
                cyano_preserved = False
                print(f"Cyano group not preserved in reaction (retrosynthetic view): {rsmi}")

            if product_has_methoxy and not reactants_have_methoxy:
                methoxy_preserved = False
                print(f"Methoxy group not preserved in reaction (retrosynthetic view): {rsmi}")

            # Record if we've found these groups
            if reactants_have_cyano or product_has_cyano:
                found_cyano = True

            if reactants_have_methoxy or product_has_methoxy:
                found_methoxy = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    def is_methoxy_present(smiles):
        """Helper function to check for methoxy groups in a molecule"""
        if not checker.check_fg("Ether", smiles):
            return False

        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Check for C-O-CH3 pattern (methoxy group)
            return mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#8]-[#6H3]"))
        return False

    # Start traversal
    dfs_traverse(route)

    # Check if functional groups are preserved throughout synthesis
    # and that they were actually present
    result = cyano_preserved and methoxy_preserved and found_cyano and found_methoxy
    print(f"Cyano preserved: {cyano_preserved}, Methoxy preserved: {methoxy_preserved}")
    print(f"Found cyano: {found_cyano}, Found methoxy: {found_methoxy}")
    print(f"Overall result: {result}")

    return result
