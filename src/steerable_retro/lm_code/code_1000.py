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
    Detects if the synthetic route uses late-stage amide formation (depth 0-1).
    """
    late_stage_amide_formation = False
    acid_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")
    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                products_mol = Chem.MolFromSmiles(products_smiles)

                if reactants_mol and products_mol and depth <= 1:
                    # Check for acid chloride in reactants and amide in products
                    if reactants_mol.HasSubstructMatch(
                        acid_chloride_pattern
                    ) and products_mol.HasSubstructMatch(amide_pattern):
                        late_stage_amide_formation = True
                        print(f"Found late-stage amide formation at depth {depth}: {rsmi}")
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late-stage amide formation detected: {late_stage_amide_formation}")
    return late_stage_amide_formation
