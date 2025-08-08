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
    Detects if the synthesis route involves protection of carboxylic acid with tert-butyl group.
    """
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a protection reaction
                if checker.check_reaction("Protection of carboxylic acid", rsmi):
                    print(f"Found potential carboxylic acid protection reaction: {rsmi}")

                    # Check if reactant has carboxylic acid
                    has_carboxylic_acid = False
                    for reactant in reactants:
                        if reactant.strip() and checker.check_fg("Carboxylic acid", reactant):
                            has_carboxylic_acid = True
                            print(f"Found carboxylic acid in reactant: {reactant}")
                            break

                    # Check if product has tert-butyl ester
                    if has_carboxylic_acid and product.strip():
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and checker.check_fg("Ester", product):
                            # Check specifically for tert-butyl ester
                            tert_butyl_ester = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                            if product_mol.HasSubstructMatch(tert_butyl_ester):
                                print(f"Found tert-butyl ester in product: {product}")
                                protection_found = True

                # Check for any reaction that converts carboxylic acid to tert-butyl ester
                if not protection_found:
                    for reactant in reactants:
                        if reactant.strip() and checker.check_fg("Carboxylic acid", reactant):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and checker.check_fg("Ester", product):
                                # Check specifically for tert-butyl ester
                                tert_butyl_ester = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                                if product_mol.HasSubstructMatch(tert_butyl_ester):
                                    print(
                                        f"Found carboxylic acid to tert-butyl ester conversion: {rsmi}"
                                    )
                                    protection_found = True
                                    break

                # Also check for Boc protection which uses tert-butyl groups
                if not protection_found and checker.check_reaction("Boc amine protection", rsmi):
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check for tert-butyloxycarbonyl (Boc) group
                        boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
                        if product_mol.HasSubstructMatch(boc_pattern):
                            print(f"Found Boc protection (tert-butyl) in product: {product}")
                            protection_found = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found
