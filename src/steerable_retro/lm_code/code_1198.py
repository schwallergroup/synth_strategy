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
    Detects if the synthesis contains a dichlorophenyl motif throughout.
    """
    has_dichlorophenyl = False

    def dfs_traverse(node):
        nonlocal has_dichlorophenyl

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for dichlorophenyl motif using SMARTS patterns for all possible isomers
                    # ortho-dichlorophenyl
                    ortho_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)cccc1")
                    # meta-dichlorophenyl
                    meta_pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)ccc1")
                    # para-dichlorophenyl
                    para_pattern = Chem.MolFromSmarts("c1c(Cl)ccc(Cl)c1")

                    # Check if the molecule contains any of the dichlorophenyl patterns
                    if (
                        mol.HasSubstructMatch(ortho_pattern)
                        or mol.HasSubstructMatch(meta_pattern)
                        or mol.HasSubstructMatch(para_pattern)
                    ):
                        print(f"Detected dichlorophenyl motif in molecule: {smiles}")
                        has_dichlorophenyl = True

                    # Also check using the checker function for aromatic halides
                    # This is a more general check that might catch other patterns
                    if checker.check_fg("Aromatic halide", smiles):
                        # Count chlorines attached to aromatic rings
                        patt = Chem.MolFromSmarts("c-Cl")
                        matches = mol.GetSubstructMatches(patt)
                        if len(matches) >= 2:
                            print(f"Detected molecule with at least 2 aromatic chlorines: {smiles}")
                            has_dichlorophenyl = True

        elif node["type"] == "reaction":
            # Check if the reaction produces a dichlorophenyl motif
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if the product contains a dichlorophenyl motif
                if product:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        # Check for dichlorophenyl motif in the product
                        ortho_pattern = Chem.MolFromSmarts("c1c(Cl)c(Cl)cccc1")
                        meta_pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)ccc1")
                        para_pattern = Chem.MolFromSmarts("c1c(Cl)ccc(Cl)c1")

                        if (
                            mol.HasSubstructMatch(ortho_pattern)
                            or mol.HasSubstructMatch(meta_pattern)
                            or mol.HasSubstructMatch(para_pattern)
                        ):
                            print(f"Detected dichlorophenyl motif in reaction product: {product}")
                            has_dichlorophenyl = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    print("Starting traversal to find dichlorophenyl motif")
    dfs_traverse(route)
    print(f"Dichlorophenyl motif found: {has_dichlorophenyl}")
    return has_dichlorophenyl
