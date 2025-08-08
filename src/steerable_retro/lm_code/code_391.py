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
    Detects if a methoxy group is preserved throughout the entire synthesis route.
    """
    # Track if we've found a methoxy group that's preserved
    target_has_methoxy = False
    all_methoxy_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal target_has_methoxy, all_methoxy_preserved

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule has a methoxy group
            # Using the checker function to detect ether and then confirming it's a methoxy
            has_methoxy = False
            if checker.check_fg("Ether", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Look for -OCH3 pattern (methoxy group)
                    methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)

            # If this is the target molecule (depth 0) and it has a methoxy group
            if depth == 0:
                if has_methoxy:
                    print(f"Target molecule has methoxy group: {mol_smiles}")
                    target_has_methoxy = True
                else:
                    print(f"Target molecule does not have methoxy group: {mol_smiles}")
                    # If target doesn't have methoxy, no need to check preservation
                    return

            # If this is a starting material, check if it has methoxy
            if node.get("in_stock", False) and has_methoxy:
                print(f"Starting material has methoxy group: {mol_smiles}")

        elif node["type"] == "reaction" and target_has_methoxy:
            # Check if the methoxy group is preserved in this reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if product has methoxy
                product_has_methoxy = False
                if checker.check_fg("Ether", product):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                        product_has_methoxy = product_mol.HasSubstructMatch(methoxy_pattern)

                # Check if at least one reactant has methoxy
                reactant_has_methoxy = False
                for r in reactants:
                    if checker.check_fg("Ether", r):
                        reactant_mol = Chem.MolFromSmiles(r)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
                        ):
                            reactant_has_methoxy = True
                            break

                # In retrosynthesis, product is what we're making and reactants are what we're breaking down to
                # If product has methoxy but none of the reactants do, methoxy is not preserved
                if product_has_methoxy and not reactant_has_methoxy:
                    print(f"Methoxy group not preserved in reaction: {rsmi}")
                    all_methoxy_preserved = False
            except Exception as e:
                print(f"Error processing reaction: {e}")
                all_methoxy_preserved = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Only return True if target has methoxy and all methoxy groups are preserved
    return target_has_methoxy and all_methoxy_preserved
