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
    Detects if the synthesis involves a sequence of functional group interconversions
    (specifically phenol→thiol, nitro→amine in this case).
    """
    phenol_to_thiol = False
    nitro_to_amine = False

    def dfs_traverse(node):
        nonlocal phenol_to_thiol, nitro_to_amine

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol to thiol conversion
            phenol_pattern = Chem.MolFromSmarts("[OH][c]")
            thiol_pattern = Chem.MolFromSmarts("[SH][c]")

            # Check for nitro to amine reduction
            nitro_pattern = Chem.MolFromSmarts("[#7+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            product_mol = Chem.MolFromSmiles(product_smiles)

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)

                if reactant_mol and product_mol:
                    # Check phenol to thiol
                    if (
                        reactant_mol.HasSubstructMatch(phenol_pattern)
                        and product_mol.HasSubstructMatch(thiol_pattern)
                        and not product_mol.HasSubstructMatch(phenol_pattern)
                    ):
                        phenol_to_thiol = True

                    # Check nitro to amine
                    if (
                        reactant_mol.HasSubstructMatch(nitro_pattern)
                        and product_mol.HasSubstructMatch(amine_pattern)
                        and not product_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        nitro_to_amine = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if both transformations are found
    return phenol_to_thiol and nitro_to_amine
