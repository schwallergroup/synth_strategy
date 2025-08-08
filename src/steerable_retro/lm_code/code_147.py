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
    This function detects if the synthetic route involves a tert-butyl ester deprotection.
    Looks for reactions where a tert-butyl ester is converted to a carboxylic acid.
    """
    found_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}: Examining reaction: {rsmi}")

            # Check if this is a deprotection reaction
            # First, check if any reactant has a tert-butyl ester group
            has_tert_butyl_ester = False
            for reactant in reactants:
                if not reactant:
                    continue

                # Check for ester functional group in reactant
                if checker.check_fg("Ester", reactant):
                    # Convert to molecule to check for tert-butyl group
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Check for tert-butyl pattern in the ester
                        tert_butyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC=O")
                        if reactant_mol.HasSubstructMatch(tert_butyl_pattern):
                            print(f"Depth {depth}: Found tert-butyl ester in reactant: {reactant}")
                            has_tert_butyl_ester = True
                            break

            # If reactant has tert-butyl ester, check if product has carboxylic acid
            if has_tert_butyl_ester and checker.check_fg("Carboxylic acid", product):
                print(f"Depth {depth}: Found carboxylic acid in product: {product}")
                # Verify this is a deprotection reaction (loss of tert-butyl group)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check that the product doesn't have the tert-butyl ester anymore
                    tert_butyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC=O")
                    if not product_mol.HasSubstructMatch(tert_butyl_pattern):
                        print(f"Depth {depth}: Confirmed tert-butyl ester deprotection")
                        found_deprotection = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_deprotection
