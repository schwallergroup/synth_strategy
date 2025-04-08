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
    Detects a convergent synthesis strategy involving amide coupling of two fragments,
    where one fragment contains a trifluoromethyl group and the other contains a bromo-aromatic ring.
    """
    # Initialize flags
    has_amide_coupling = False
    has_trifluoromethyl = False
    has_bromo_aromatic = False

    def dfs_traverse(node):
        nonlocal has_amide_coupling, has_trifluoromethyl, has_bromo_aromatic

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide coupling (C-N bond formation)
            if len(reactants) >= 2:  # Multiple reactants indicate potential convergent step
                # Check if one reactant has an acid chloride (C(=O)Cl) and another has an amine
                has_acid_chloride = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[Cl]")):
                            has_acid_chloride = True
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                            has_amine = True

                # Check if product has an amide bond
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[NH]")
                ):
                    if has_acid_chloride and has_amine:
                        print("Found amide coupling reaction")
                        has_amide_coupling = True

            # Check for trifluoromethyl group in reactants or product
            for smi in reactants + [product]:
                mol = Chem.MolFromSmiles(smi)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]([F])([F])[F]")):
                    print("Found trifluoromethyl group")
                    has_trifluoromethyl = True

                # Check for bromo-aromatic in reactants or product
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br]")):
                    print("Found bromo-aromatic group")
                    has_bromo_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if all required elements are found
    return has_amide_coupling and has_trifluoromethyl and has_bromo_aromatic
