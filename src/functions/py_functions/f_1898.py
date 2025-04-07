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
    This function detects late-stage amide coupling (depth 0-1) connecting complex fragments.
    """
    late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0-1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation
                if len(reactants_smiles) >= 2:  # At least two reactants
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[#8])")
                    amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                    amide_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#7]")

                    # Check reactants for acid and amine
                    has_acid = False
                    has_amine = False

                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                has_acid = True
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    # Check product for amide
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    has_amide = product_mol and product_mol.HasSubstructMatch(
                        amide_pattern
                    )

                    if has_acid and has_amine and has_amide:
                        print(f"Late-stage amide coupling detected at depth {depth}")
                        late_stage_amide_coupling = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_coupling
