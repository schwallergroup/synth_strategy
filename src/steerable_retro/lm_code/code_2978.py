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
    Detects if the synthesis route involves a late-stage C-C bond formation
    between two aromatic fragments.
    """
    late_stage_cc_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cc_formation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Parse molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol and len(reactant_mols) >= 2:
                    # Check if both reactants have aromatic rings
                    aromatic_reactants = 0
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1")):
                            aromatic_reactants += 1

                    # If we have at least 2 aromatic reactants and the product has a new C-C bond
                    # between aromatic rings, this is likely a late-stage aromatic C-C bond formation
                    if aromatic_reactants >= 2:
                        print(f"Found late-stage aromatic C-C bond formation at depth {depth}")
                        late_stage_cc_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_cc_formation
