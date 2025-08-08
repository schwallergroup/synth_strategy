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
    This function detects the use of tert-butyl protection for carboxylic acids
    in the synthetic route.
    """
    protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a protection or deprotection reaction directly
            if checker.check_reaction("Protection of carboxylic acid", rsmi):
                print(f"Found carboxylic acid protection reaction: {rsmi}")
                protection_found = True
                return

            if checker.check_reaction("Deprotection of carboxylic acid", rsmi):
                print(f"Found carboxylic acid deprotection reaction: {rsmi}")
                protection_found = True
                return

            if checker.check_reaction("COOH ethyl deprotection", rsmi):
                print(f"Found COOH deprotection reaction: {rsmi}")
                # Need to verify it's specifically a tert-butyl deprotection
                for reactant in reactants:
                    if not reactant:
                        continue
                    try:
                        if checker.check_fg("Ester", reactant):
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol:
                                tert_butyl_pattern = Chem.MolFromSmarts(
                                    "[O;D2][C](=[O])[C]-[C]([C])([C])[C]"
                                )
                                if react_mol.HasSubstructMatch(tert_butyl_pattern):
                                    print(f"Found tert-butyl ester in reactant: {reactant}")
                                    protection_found = True
                                    return
                    except Exception as e:
                        print(f"Error checking reactant for tert-butyl ester: {e}")

            # Forward direction: Carboxylic acid in reactants, tert-butyl ester in product
            carboxylic_acid_in_reactants = False
            for reactant in reactants:
                if reactant and checker.check_fg("Carboxylic acid", reactant):
                    print(f"Found carboxylic acid in reactant: {reactant}")
                    carboxylic_acid_in_reactants = True
                    break

            if carboxylic_acid_in_reactants and product:
                try:
                    # Check for tert-butyl ester in product
                    if checker.check_fg("Ester", product):
                        # Verify it's specifically a tert-butyl ester
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            tert_butyl_pattern = Chem.MolFromSmarts(
                                "[O;D2][C](=[O])[C]-[C]([C])([C])[C]"
                            )
                            if prod_mol.HasSubstructMatch(tert_butyl_pattern):
                                print(f"Found tert-butyl ester in product: {product}")
                                protection_found = True
                                return
                except Exception as e:
                    print(f"Error checking product for tert-butyl ester: {e}")

            # Retrosynthetic direction: tert-butyl ester in reactants, carboxylic acid in product
            tert_butyl_in_reactants = False
            for reactant in reactants:
                if not reactant:
                    continue
                try:
                    if checker.check_fg("Ester", reactant):
                        # Verify it's specifically a tert-butyl ester
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            tert_butyl_pattern = Chem.MolFromSmarts(
                                "[O;D2][C](=[O])[C]-[C]([C])([C])[C]"
                            )
                            if react_mol.HasSubstructMatch(tert_butyl_pattern):
                                print(f"Found tert-butyl ester in reactant: {reactant}")
                                tert_butyl_in_reactants = True
                                break
                except Exception as e:
                    print(f"Error checking reactant for tert-butyl ester: {e}")

            if tert_butyl_in_reactants and product and checker.check_fg("Carboxylic acid", product):
                print(
                    f"Found carboxylic acid in product with tert-butyl ester in reactant: {product}"
                )
                protection_found = True
                return

            # Direct check for tert-butyl groups in the reaction
            for reactant in reactants:
                if not reactant:
                    continue
                react_mol = Chem.MolFromSmiles(reactant)
                if react_mol:
                    # Check for tert-butyl group connected to oxygen (tert-butyl ester or other)
                    tert_butyl_pattern = Chem.MolFromSmarts("[O]-[C]([C])([C])[C]")
                    if react_mol.HasSubstructMatch(tert_butyl_pattern):
                        print(f"Found tert-butyl group in reactant: {reactant}")
                        # Check if this is part of an ester
                        if checker.check_fg("Ester", reactant):
                            print(f"Confirmed tert-butyl ester in reactant: {reactant}")
                            protection_found = True
                            return

            # Check product for tert-butyl groups
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol:
                tert_butyl_pattern = Chem.MolFromSmarts("[O]-[C]([C])([C])[C]")
                if prod_mol.HasSubstructMatch(tert_butyl_pattern):
                    print(f"Found tert-butyl group in product: {product}")
                    # Check if this is part of an ester
                    if checker.check_fg("Ester", product):
                        print(f"Confirmed tert-butyl ester in product: {product}")
                        protection_found = True
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if not protection_found:
        print("No tert-butyl protection of carboxylic acid found in the route")

    return protection_found
