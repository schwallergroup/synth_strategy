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


def main(route, depth=0, max_depth=2):
    """
    Detects if the synthesis route includes a late-stage N-oxidation.
    Late-stage means it occurs within the first few steps (low depth in retrosynthetic analysis).
    """
    if route["type"] == "reaction" and depth <= max_depth:
        # Check if this is an N-oxidation reaction
        try:
            rsmi = route["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains N-oxide but reactants don't
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                has_n_oxide_product = False
                for atom in product_mol.GetAtoms():
                    if (
                        atom.GetSymbol() == "N"
                        and atom.GetFormalCharge() == 1
                        and any(n.GetSymbol() == "O" for n in atom.GetNeighbors())
                    ):
                        has_n_oxide_product = True
                        break

                if has_n_oxide_product:
                    # Check if reactants don't have N-oxide
                    has_n_oxide_reactant = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            for atom in reactant_mol.GetAtoms():
                                if (
                                    atom.GetSymbol() == "N"
                                    and atom.GetFormalCharge() == 1
                                    and any(
                                        n.GetSymbol() == "O"
                                        for n in atom.GetNeighbors()
                                    )
                                ):
                                    has_n_oxide_reactant = True
                                    break

                    if not has_n_oxide_reactant:
                        print(f"Late-stage N-oxidation detected at depth {depth}")
                        return True
        except Exception as e:
            print(f"Error analyzing N-oxidation: {e}")

    # Recursively check children
    for child in route.get("children", []):
        if late_stage_n_oxidation_strategy(child, depth + 1, max_depth):
            return True

    return False
