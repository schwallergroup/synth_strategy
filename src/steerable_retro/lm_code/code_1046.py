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
    Detects if the synthetic route involves sequential halogen manipulations,
    such as halogen exchange or replacement with other functional groups.
    """
    halogen_transformations = 0

    def dfs_traverse(node):
        nonlocal halogen_transformations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for halogen exchange or replacement
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Check for aryl halide in reactant
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,Cl,I]")):
                            # Case 1: Halogen exchange (one halogen to another)
                            if product_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[Br,Cl,I]")
                            ) and not reactant_mol.GetSubstructMatch(
                                Chem.MolFromSmarts("[c]-[Br]")
                            ) == product_mol.GetSubstructMatch(
                                Chem.MolFromSmarts("[c]-[Br]")
                            ):
                                print("Found halogen exchange")
                                halogen_transformations += 1

                            # Case 2: Halogen replacement with other group (O, N, B, etc.)
                            elif not product_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[Br,Cl,I]")
                            ) and (
                                product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[O,N,B,S]"))
                                or product_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[c]"))
                            ):
                                print("Found halogen replacement")
                                halogen_transformations += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogen_transformations >= 2  # At least two halogen transformations
