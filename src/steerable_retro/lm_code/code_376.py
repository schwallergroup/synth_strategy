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
    This function detects late-stage C-N bond formations including reductive amination
    and aryl-N coupling.
    """
    # Track C-N bond formations
    reductive_amination = False
    aryl_n_coupling = False
    late_stage_cn_formation = False

    # SMARTS patterns
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
    amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(N#*)]")
    aryl_pattern = Chem.MolFromSmarts("c")

    def dfs_traverse(node):
        nonlocal reductive_amination, aryl_n_coupling, late_stage_cn_formation

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check for reductive amination (carbonyl + amine → amine)
                if (
                    reactants_mol.HasSubstructMatch(carbonyl_pattern)
                    and reactants_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                    and not product_mol.HasSubstructMatch(carbonyl_pattern)
                ):
                    reductive_amination = True
                    # Check if this is a late-stage reaction (depth 0 or 1)
                    if node["metadata"].get("depth", 0) <= 1:
                        late_stage_cn_formation = True
                        print("Detected late-stage reductive amination")

                # Check for aryl-N coupling
                if (
                    reactants_mol.HasSubstructMatch(aryl_pattern)
                    and reactants_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(aryl_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    # This is a simplistic check - in a real implementation,
                    # we would need to compare atom mappings to confirm new C-N bond
                    aryl_n_coupling = True
                    # Check if this is a late-stage reaction (depth 0 or 1)
                    if node["metadata"].get("depth", 0) <= 1:
                        late_stage_cn_formation = True
                        print("Detected late-stage aryl-N coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have late-stage C-N bond formation
    strategy_present = late_stage_cn_formation and (reductive_amination or aryl_n_coupling)

    if strategy_present:
        print("Late-stage C-N bond formation strategy detected")

    return strategy_present
