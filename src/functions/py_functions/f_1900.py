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
    This function detects early-stage oxime formation from ketone.
    """
    oxime_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal oxime_formation_detected

        if node["type"] == "reaction" and depth >= 4:  # Early stage (depth 4+)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for oxime formation
                ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")
                oxime_pattern = Chem.MolFromSmarts("[#6]=[#7]-[#8]")
                hydroxylamine_pattern = Chem.MolFromSmarts("[#7]-[#8]")

                # Check if ketone is in reactants
                ketone_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(ketone_pattern):
                        ketone_in_reactants = True
                        break

                # Check if oxime is in product and hydroxylamine is in reactants
                product_mol = Chem.MolFromSmiles(product_smiles)
                if (
                    ketone_in_reactants
                    and product_mol
                    and product_mol.HasSubstructMatch(oxime_pattern)
                ):
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            hydroxylamine_pattern
                        ):
                            print(f"Early oxime formation detected at depth {depth}")
                            oxime_formation_detected = True
                            break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return oxime_formation_detected
