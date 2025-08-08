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
    Detects ketone to carboxylic acid transformation in the synthesis route.
    """
    transformation_found = False

    def dfs_traverse(node):
        nonlocal transformation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for C(=O)CH3 to C(=O)OH transformation
            if "[C:" in product and "=[O:" in product and "[OH:" in product:
                product_mol = Chem.MolFromSmiles(product)

                # Look for carboxylic acid in product
                for atom in product_mol.GetAtoms():
                    if atom.GetSymbol() == "C":
                        oh_neighbor = None
                        double_o_neighbor = None

                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetSymbol() == "O":
                                bond = product_mol.GetBondBetweenAtoms(
                                    atom.GetIdx(), neighbor.GetIdx()
                                )
                                if (
                                    bond.GetBondType() == Chem.BondType.SINGLE
                                    and neighbor.GetTotalNumHs() > 0
                                ):
                                    oh_neighbor = neighbor
                                elif bond.GetBondType() == Chem.BondType.DOUBLE:
                                    double_o_neighbor = neighbor

                        if oh_neighbor and double_o_neighbor:
                            # Check if there was a methyl ketone in reactants
                            for r in reactants:
                                if "[C:" in r and "=[O:" in r and "[CH3:" in r:
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol:
                                        for r_atom in r_mol.GetAtoms():
                                            if r_atom.GetSymbol() == "C":
                                                methyl_neighbor = None
                                                double_o_neighbor = None

                                                for neighbor in r_atom.GetNeighbors():
                                                    if (
                                                        neighbor.GetSymbol() == "C"
                                                        and neighbor.GetTotalNumHs() == 3
                                                    ):
                                                        methyl_neighbor = neighbor
                                                    elif (
                                                        neighbor.GetSymbol() == "O"
                                                        and r_mol.GetBondBetweenAtoms(
                                                            r_atom.GetIdx(), neighbor.GetIdx()
                                                        ).GetBondType()
                                                        == Chem.BondType.DOUBLE
                                                    ):
                                                        double_o_neighbor = neighbor

                                                if methyl_neighbor and double_o_neighbor:
                                                    transformation_found = True
                                                    print(
                                                        "Ketone to carboxylic acid transformation detected"
                                                    )
                                                    return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_found
