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
    This function detects an ester hydrolysis followed by amide formation sequence.
    """
    ester_hydrolysis_reactions = []
    amide_formation_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for ester hydrolysis
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")
                acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[#8])")

                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for ester in reactants and acid in product
                if product_mol and product_mol.HasSubstructMatch(acid_pattern):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            ester_pattern
                        ):
                            ester_hydrolysis_reactions.append((depth, node))
                            print(f"Ester hydrolysis detected at depth {depth}")
                            break

                # Check for amide formation
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#7]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    has_acid = False
                    has_amine = False

                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(acid_pattern):
                                has_acid = True
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_acid and has_amine:
                        amide_formation_reactions.append((depth, node))
                        print(f"Amide formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if ester hydrolysis is followed by amide formation
    if ester_hydrolysis_reactions and amide_formation_reactions:
        for ester_depth, _ in ester_hydrolysis_reactions:
            for amide_depth, _ in amide_formation_reactions:
                if amide_depth < ester_depth:  # Lower depth means later in synthesis
                    print(
                        "Ester hydrolysis followed by amide formation strategy detected"
                    )
                    return True

    return False
