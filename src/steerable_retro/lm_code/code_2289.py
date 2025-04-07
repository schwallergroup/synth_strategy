#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects late-stage carboxylic acid protection strategy where
    a carboxylic acid is converted to a tert-butyl ester in the final steps.
    In retrosynthesis, we look for the reverse: tert-butyl ester being deprotected to carboxylic acid.
    """
    protection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
            and depth <= 1
        ):
            # Only consider reactions at depth 0 or 1 (late stage)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # In retrosynthesis, product would have carboxylic acid, reactants would have ester
            product_has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
            if product_has_carboxylic_acid:
                print(f"Found carboxylic acid in product: {product}")

            # Check if any reactant contains an ester
            reactants_with_ester = []
            for reactant in reactants:
                if checker.check_fg("Ester", reactant):
                    reactants_with_ester.append(reactant)
                    print(f"Found ester in reactant: {reactant}")

            if product_has_carboxylic_acid and reactants_with_ester:
                # Check if this is a deprotection reaction (protection in forward direction)
                deprotection_reactions = [
                    "COOH ethyl deprotection",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Deprotection of carboxylic acid",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                ]

                for reaction_type in deprotection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found deprotection reaction: {reaction_type}")

                        # Check if any reactant contains a tert-butyl ester
                        for reactant in reactants_with_ester:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Check for tert-butyl group connected to ester
                                tert_butyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                                if reactant_mol.HasSubstructMatch(tert_butyl_pattern):
                                    print(
                                        f"Late-stage carboxylic acid protection (tert-butyl ester) detected at depth {depth}"
                                    )
                                    protection_detected = True
                                    break

                                # Alternative pattern for tert-butyl ester
                                alt_pattern = Chem.MolFromSmarts("CC(C)(C)O[C,c]=O")
                                if reactant_mol.HasSubstructMatch(alt_pattern):
                                    print(
                                        f"Late-stage carboxylic acid protection (tert-butyl ester) detected at depth {depth}"
                                    )
                                    protection_detected = True
                                    break

                        if not protection_detected:
                            print("Deprotection found but not with tert-butyl group")

                if not protection_detected:
                    # Check if it's a general esterification/deprotection reaction
                    for reactant in reactants_with_ester:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Check for tert-butyl group connected to ester
                            tert_butyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                            alt_pattern = Chem.MolFromSmarts("CC(C)(C)O[C,c]=O")
                            if reactant_mol.HasSubstructMatch(
                                tert_butyl_pattern
                            ) or reactant_mol.HasSubstructMatch(alt_pattern):
                                print(
                                    f"Late-stage carboxylic acid protection (tert-butyl ester) detected at depth {depth}"
                                )
                                protection_detected = True
                                break
            else:
                if not product_has_carboxylic_acid:
                    print("No carboxylic acid found in product")
                if not reactants_with_ester:
                    print("No ester found in reactants")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Protection detected: {protection_detected}")
    return protection_detected
