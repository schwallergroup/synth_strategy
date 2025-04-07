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


def main(route):
    """
    Detects if the synthetic route includes a protection-deprotection sequence for carboxylic acid
    (esterification followed by hydrolysis).
    """
    # Track reactions by depth
    reactions_by_depth = {}

    def collect_reactions(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            if depth not in reactions_by_depth:
                reactions_by_depth[depth] = []
            reactions_by_depth[depth].append(node["metadata"]["rsmi"])

        for child in node.get("children", []):
            collect_reactions(child, depth + 1)

    collect_reactions(route)

    # Check for esterification followed by hydrolysis
    esterification_depths = []
    hydrolysis_depths = []

    for depth, reactions in reactions_by_depth.items():
        for rsmi in reactions:
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for esterification
                has_acid_reactant = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[OH]"))
                    for r in reactant_mols
                )
                has_ester_product = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[O][C]")
                )
                if has_acid_reactant and has_ester_product:
                    esterification_depths.append(depth)

                # Check for hydrolysis
                has_ester_reactant = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][C]"))
                    for r in reactant_mols
                )
                has_acid_product = product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[OH]")
                )
                if has_ester_reactant and has_acid_product:
                    hydrolysis_depths.append(depth)
            except:
                print(f"Error processing reaction SMILES at depth {depth}")

    # Check if there's an esterification at a higher depth followed by hydrolysis at a lower depth
    for ester_depth in esterification_depths:
        for hydro_depth in hydrolysis_depths:
            if (
                ester_depth > hydro_depth
            ):  # Remember higher depth = earlier in synthesis
                print(
                    f"Protection-deprotection sequence detected: esterification at depth {ester_depth}, hydrolysis at depth {hydro_depth}"
                )
                return True

    return False
