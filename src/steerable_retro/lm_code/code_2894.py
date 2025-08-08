#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects a late-stage amide coupling in the synthesis route.
    """
    late_stage_amide_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage = low depth
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants include amine and carboxylic acid derivative
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                amine_pattern = Chem.MolFromSmarts("[#6]-[N;H2]")
                carboxylic_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
                amide_pattern = Chem.MolFromSmarts("[#6]-[N;H]-[C](=[O])")

                reactant_has_amine = any(
                    mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                )
                reactant_has_carboxylic = any(
                    mol and mol.HasSubstructMatch(carboxylic_pattern) for mol in reactant_mols
                )
                product_has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if reactant_has_amine and reactant_has_carboxylic and product_has_amide:
                    print("Late-stage amide coupling detected")
                    late_stage_amide_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_coupling_found
