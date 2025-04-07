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
    Detects convergent synthesis with late-stage etherification connecting two complex fragments.
    """
    late_stage_etherification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_etherification

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Late-stage reaction (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (convergent)
                if len(reactants) >= 2:
                    # Check for etherification (C-O-C formation)
                    prod_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Look for ether pattern in product
                    ether_patt = Chem.MolFromSmarts("[#6][#8][#6]")
                    if prod_mol and prod_mol.HasSubstructMatch(ether_patt):
                        # Check if ether bond is newly formed
                        ether_count_prod = len(prod_mol.GetSubstructMatches(ether_patt))
                        ether_count_reactants = sum(
                            len(r.GetSubstructMatches(ether_patt)) if r else 0
                            for r in reactant_mols
                        )

                        if ether_count_prod > ether_count_reactants:
                            print(f"Found late-stage etherification at depth {depth}")
                            late_stage_etherification = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_etherification
