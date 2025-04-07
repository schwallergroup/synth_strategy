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
    Detects if the synthesis uses a late-stage amide coupling strategy
    (amide formation in the final or penultimate step)
    """
    amide_coupling_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_at_late_stage

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction
                # Look for C(=O)-N bond formation
                prod_mol = Chem.MolFromSmiles(product)

                # Check if any reactant has an amine group
                has_amine = False
                has_carboxyl_derivative = False

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Check for amine
                        if r_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2,NH;!$(NC=O)]")):
                            has_amine = True
                        # Check for carboxylic acid derivative
                        if r_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6](=[#8])[#8,#7,Cl,F,Br,I]")
                        ):
                            has_carboxyl_derivative = True

                # Check if product has an amide bond
                if prod_mol and prod_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6](=[#8])[#7;!$(N=*)]")
                ):
                    if has_amine and has_carboxyl_derivative:
                        print(f"Found amide coupling at depth {depth}")
                        amide_coupling_at_late_stage = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return amide_coupling_at_late_stage
