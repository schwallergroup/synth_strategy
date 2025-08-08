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
    This function detects if the synthetic route employs a Suzuki coupling strategy.
    It looks for reactions where a bromo-aromatic and boronic acid derivative form a biaryl C-C bond.
    """
    found_suzuki = False

    def dfs_traverse(node):
        nonlocal found_suzuki

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains bromo-aromatic
            bromo_pattern = Chem.MolFromSmarts("[c]-[Br]")
            # Check if one reactant contains boronic acid derivative
            boronic_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
            # Check if product contains biaryl bond that wasn't in reactants
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_bromo = False
            has_boronic = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(bromo_pattern):
                        has_bromo = True
                    if mol and mol.HasSubstructMatch(boronic_pattern):
                        has_boronic = True
                except:
                    continue

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if (
                    has_bromo
                    and has_boronic
                    and prod_mol
                    and prod_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    print("Found Suzuki coupling reaction:", rsmi)
                    found_suzuki = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_suzuki
