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
    Detects a synthetic strategy involving building a carbon framework first,
    followed by functional group transformations (reduction, leaving group installation).
    """
    # Track reaction types and their depths
    c_c_bond_formation_depth = -1
    reduction_depth = -1
    leaving_group_installation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal c_c_bond_formation_depth, reduction_depth, leaving_group_installation_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if reactants and product:
                    # Check for C-C bond formation (malonate alkylation)
                    malonate_pattern = Chem.MolFromSmarts("[CH2][C](=[O])[O][CH2][CH3]")
                    alkyl_halide_pattern = Chem.MolFromSmarts("[Br,Cl,I][CH2]")
                    alkylated_malonate_pattern = Chem.MolFromSmarts(
                        "[CH]([CH2][*])([C](=[O])[O][CH2][CH3])[C](=[O])[O][CH2][CH3]"
                    )

                    has_malonate = False
                    has_alkyl_halide = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(malonate_pattern):
                                has_malonate = True
                            if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        if product_mol.HasSubstructMatch(alkylated_malonate_pattern):
                            if has_malonate and has_alkyl_halide:
                                c_c_bond_formation_depth = depth

                        # Check for reduction (diester to diol)
                        diester_pattern = Chem.MolFromSmarts("[C](=[O])[O][CH2][CH3]")
                        diol_pattern = Chem.MolFromSmarts("[OH][CH2][CH][CH2][OH]")

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(diester_pattern):
                                if product_mol.HasSubstructMatch(diol_pattern):
                                    reduction_depth = depth

                        # Check for leaving group installation (alcohol to tosylate)
                        alcohol_pattern = Chem.MolFromSmarts("[OH][CH2]")
                        tosyl_chloride_pattern = Chem.MolFromSmarts("Cl[S](=[O])(=[O])[c]")
                        tosylate_pattern = Chem.MolFromSmarts("[O][S](=[O])(=[O])[c]")

                        has_alcohol = False
                        has_tosyl_chloride = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(alcohol_pattern):
                                    has_alcohol = True
                                if reactant_mol.HasSubstructMatch(tosyl_chloride_pattern):
                                    has_tosyl_chloride = True

                        if product_mol.HasSubstructMatch(tosylate_pattern):
                            if has_alcohol and has_tosyl_chloride:
                                leaving_group_installation_depth = depth

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy follows the pattern: C-C bond formation -> reduction -> leaving group installation
    # Note: Lower depth values are later in the synthesis (closer to final product)
    has_build_and_functionalize = (
        c_c_bond_formation_depth > reduction_depth > leaving_group_installation_depth >= 0
    )

    print(f"Build and functionalize strategy: {has_build_and_functionalize}")
    print(f"  C-C bond formation depth: {c_c_bond_formation_depth}")
    print(f"  Reduction depth: {reduction_depth}")
    print(f"  Leaving group installation depth: {leaving_group_installation_depth}")

    return has_build_and_functionalize
