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
    This function detects if the synthesis involves construction of a heterocycle
    through imine cyclization.
    """
    has_imine_cyclization = False

    def dfs_traverse(node):
        nonlocal has_imine_cyclization

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketone and amine/sulfonamide in reactants
            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C,N]")
            amine_pattern = Chem.MolFromSmarts("[N][H]")
            sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=[O])(=[O])")

            # Check for cyclic imine in product
            cyclic_imine_pattern = Chem.MolFromSmarts("[C]=[N]")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(cyclic_imine_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(ketone_pattern) and (
                            reactant_mol.HasSubstructMatch(amine_pattern)
                            or reactant_mol.HasSubstructMatch(sulfonamide_pattern)
                        ):

                            # Check if product has more rings than reactants (cyclization)
                            product_rings = Chem.GetSSSR(product_mol)
                            reactant_rings = Chem.GetSSSR(reactant_mol)

                            if len(product_rings) > len(reactant_rings):
                                has_imine_cyclization = True
                                print("Detected heterocycle formation via imine cyclization")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_imine_cyclization
