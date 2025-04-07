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
    Detects incorporation of an amino acid derivative with preserved stereochemistry.
    """
    amino_acid_found = False
    stereocenter_preserved = False

    def dfs_traverse(node, depth=0):
        nonlocal amino_acid_found, stereocenter_preserved

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amino acid pattern
                amino_acid_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#6](-[*])-[#8]")

                # Check for stereochemistry in reactants
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(amino_acid_pattern):
                        amino_acid_found = True
                        # Count stereocenters in reactant
                        r_chiral_centers = len(
                            Chem.FindMolChiralCenters(r_mol, includeUnassigned=False)
                        )

                        # Check if product preserves stereocenter
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol:
                            p_chiral_centers = len(
                                Chem.FindMolChiralCenters(p_mol, includeUnassigned=False)
                            )
                            if p_chiral_centers >= r_chiral_centers and r_chiral_centers > 0:
                                stereocenter_preserved = True
                                print(
                                    f"Found amino acid with preserved stereochemistry at depth {depth}"
                                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return amino_acid_found and stereocenter_preserved
