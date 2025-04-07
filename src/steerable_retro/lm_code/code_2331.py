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
    This function detects O-arylation as a fragment coupling strategy.
    Looks for Ar-Cl + Ar-OH → Ar-O-Ar pattern.
    """
    o_arylation_found = False

    def dfs_traverse(node):
        nonlocal o_arylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have exactly 2 reactants (convergent step)
                if len(reactants) == 2:
                    # Check for aryl chloride and phenol patterns
                    has_aryl_chloride = False
                    has_phenol = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("c-Cl")):
                                has_aryl_chloride = True
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("c-[OH]")):
                                has_phenol = True

                    # Check if product has diaryl ether pattern
                    prod_mol = Chem.MolFromSmiles(product)
                    if has_aryl_chloride and has_phenol and prod_mol:
                        if prod_mol.HasSubstructMatch(Chem.MolFromSmarts("c-[O]-c")):
                            o_arylation_found = True
                            print("O-arylation fragment coupling detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return o_arylation_found
