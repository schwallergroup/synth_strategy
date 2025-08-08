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
    This function detects if the synthetic route involves coupling of a fluorinated aryl group.
    Looks for cross-coupling reactions where one reactant contains fluorine atoms on an aromatic ring.
    """
    found_fluorinated_coupling = False

    def dfs_traverse(node):
        nonlocal found_fluorinated_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains a fluorinated aryl group
            has_fluorinated_aryl = False
            has_coupling_partner = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for fluorinated aryl pattern
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][F]")):
                            has_fluorinated_aryl = True
                        # Check for coupling partner (halide or metal)
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c][Br,I,Cl]")
                        ) or mol.HasSubstructMatch(Chem.MolFromSmarts("[c][B,Sn,Zn]")):
                            has_coupling_partner = True

            # Check if product contains a biaryl system with fluorine
            if has_fluorinated_aryl and has_coupling_partner:
                prod_mol = Chem.MolFromSmiles(product)
                if (
                    prod_mol
                    and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[c]"))
                    and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][F]"))
                ):
                    print("Found fluorinated aryl coupling")
                    found_fluorinated_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_fluorinated_coupling
