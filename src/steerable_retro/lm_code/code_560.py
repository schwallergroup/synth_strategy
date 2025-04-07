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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects whether the core scaffold is preserved throughout the synthesis,
    with only functional group modifications occurring.
    """
    # Track all molecule nodes with their scaffolds and depths
    molecule_scaffolds = []

    def extract_scaffold(smiles):
        """Extract the Murcko scaffold from a molecule"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Could not parse SMILES: {smiles}")
            return None

        # Generate a Murcko scaffold (focusing on rings and connectivity)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold.GetNumAtoms() == 0:
            print(f"No scaffold found for: {smiles}")
            return None

        scaffold_smiles = Chem.MolToSmiles(scaffold)
        print(f"Extracted scaffold: {scaffold_smiles} from {smiles}")
        return scaffold

    def dfs_traverse(node, depth=0):
        """Traverse the synthesis route and collect scaffolds"""
        if node["type"] == "mol":
            # For molecule nodes, extract and store the scaffold
            scaffold = extract_scaffold(node["smiles"])
            if scaffold:
                molecule_scaffolds.append((depth, scaffold, node["smiles"]))

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort scaffolds by depth (late-stage to early-stage)
    molecule_scaffolds.sort(key=lambda x: x[0])

    # Check if we have enough molecules to compare
    if len(molecule_scaffolds) < 2:
        print("Not enough molecules to compare scaffolds")
        return False

    # The target molecule (product) scaffold is our reference
    reference_scaffold = molecule_scaffolds[0][1]
    reference_smiles = molecule_scaffolds[0][2]
    print(f"Reference scaffold from: {reference_smiles}")

    # Compare each scaffold to the reference
    scaffold_preserved = True
    for depth, scaffold, smiles in molecule_scaffolds[1:]:
        # Compare fingerprints of scaffolds
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(reference_scaffold, 2)
        curr_fp = AllChem.GetMorganFingerprintAsBitVect(scaffold, 2)

        # Calculate Tanimoto similarity
        similarity = AllChem.DataStructs.TanimotoSimilarity(ref_fp, curr_fp)
        print(f"Scaffold similarity at depth {depth}: {similarity:.2f} for {smiles}")

        # Check if the scaffolds are significantly different
        if similarity < 0.8:  # Increased threshold for scaffold preservation
            print(f"Scaffold not preserved at depth {depth}, similarity: {similarity:.2f}")
            scaffold_preserved = False
            break

    # Check if reactions only modify functional groups
    def check_reactions(node):
        """Check if reactions only modify functional groups"""
        nonlocal scaffold_preserved

        if node["type"] == "reaction" and scaffold_preserved:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a functional group transformation
                # by comparing scaffolds of reactants and products
                product_scaffold = extract_scaffold(product)
                if not product_scaffold:
                    return

                # Check if the reaction is a known functional group transformation
                is_fg_transformation = False

                # List of common functional group transformation reactions
                fg_transformation_reactions = [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Reduction of aldehydes and ketones to alcohols",
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Alcohol protection with silyl ethers",
                    "Alcohol deprotection from silyl ethers",
                    "Boc amine protection",
                    "Boc amine deprotection",
                ]

                for rxn_type in fg_transformation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected functional group transformation: {rxn_type}")
                        is_fg_transformation = True
                        break

                # If not a known FG transformation, check scaffold similarity
                if not is_fg_transformation:
                    for reactant in reactants:
                        reactant_scaffold = extract_scaffold(reactant)
                        if not reactant_scaffold:
                            continue

                        # Compare scaffolds
                        prod_fp = AllChem.GetMorganFingerprintAsBitVect(product_scaffold, 2)
                        react_fp = AllChem.GetMorganFingerprintAsBitVect(reactant_scaffold, 2)
                        similarity = AllChem.DataStructs.TanimotoSimilarity(prod_fp, react_fp)

                        print(f"Reaction scaffold similarity: {similarity:.2f}")

                        # If any reactant has a significantly different scaffold,
                        # this is not just a functional group modification
                        if similarity < 0.8:
                            print(f"Reaction modifies scaffold, similarity: {similarity:.2f}")
                            print(f"Reactant: {reactant}")
                            print(f"Product: {product}")

                            # Check if functional groups changed
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            product_mol = Chem.MolFromSmiles(product)

                            # List of functional groups to check
                            fg_list = [
                                "Aldehyde",
                                "Ketone",
                                "Carboxylic acid",
                                "Ester",
                                "Alcohol",
                                "Amine",
                                "Amide",
                                "Nitrile",
                                "Nitro group",
                                "Halide",
                                "Ether",
                                "Phenol",
                            ]

                            # Check if functional groups changed between reactant and product
                            fg_changed = False
                            for fg in fg_list:
                                reactant_has_fg = checker.check_fg(fg, reactant)
                                product_has_fg = checker.check_fg(fg, product)

                                if reactant_has_fg != product_has_fg:
                                    print(f"Functional group change detected: {fg}")
                                    fg_changed = True

                            # If no functional group changes detected, this is not a FG transformation
                            if not fg_changed:
                                scaffold_preserved = False
                                return

        # Continue traversal
        for child in node.get("children", []):
            check_reactions(child)

    # Check reactions if scaffolds are preserved
    if scaffold_preserved:
        check_reactions(route)

    print(f"Scaffold preservation strategy result: {scaffold_preserved}")
    return scaffold_preserved
