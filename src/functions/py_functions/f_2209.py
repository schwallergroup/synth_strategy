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
    Detects convergent synthesis with late-stage amide coupling.
    A convergent synthesis has multiple branches that are combined in a late stage.
    """
    is_convergent = False
    has_late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, has_late_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Check if this is a convergent step (has multiple reactants)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling
                if len(reactants) >= 2:
                    # Convert to RDKit molecules
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                        # Check for amide formation
                        amide_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[N:3]")
                        product_matches = product_mol.GetSubstructMatches(amide_pattern)

                        # Check if amide bond is newly formed
                        for match in product_matches:
                            c_idx, o_idx, n_idx = match
                            amide_bond_found = False

                            for reactant in reactant_mols:
                                if reactant.HasSubstructMatch(
                                    Chem.MolFromSmarts(
                                        f"[C:{c_idx+1}](=[O:{o_idx+1}])[N:{n_idx+1}]"
                                    )
                                ):
                                    amide_bond_found = True
                                    break

                            if not amide_bond_found:
                                has_late_amide_coupling = True
                                is_convergent = True
                                print(
                                    f"Found late-stage amide coupling at depth {depth}"
                                )
                                break
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_convergent and has_late_amide_coupling
