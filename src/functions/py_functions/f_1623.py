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
    Detects early-stage coupling of two fragments via propargyl ether formation.
    Looks for C-O bond formation between an aromatic ring and a propargyl group.
    """
    found_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling

        if node["type"] == "reaction" and depth >= 2:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a propargyl ether formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    propargyl_ether_pattern = Chem.MolFromSmarts(
                        "[c]-[O]-[CH2]-[C]#[C]-[c]"
                    )
                    if product_mol.HasSubstructMatch(propargyl_ether_pattern):
                        # Check if reactants include phenol and propargyl derivative
                        phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")
                        propargyl_pattern = Chem.MolFromSmarts(
                            "[Br,Cl,I]-[CH2]-[C]#[C]"
                        )

                        has_phenol = False
                        has_propargyl = False

                        for reactant in reactants:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol:
                                if r_mol.HasSubstructMatch(phenol_pattern):
                                    has_phenol = True
                                if r_mol.HasSubstructMatch(propargyl_pattern):
                                    has_propargyl = True

                        if has_phenol and has_propargyl:
                            found_coupling = True
                            print(
                                f"Found early-stage propargyl ether coupling at depth {depth}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_coupling
