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
    This function detects if the synthesis follows the specific sequence:
    Reduction → Disconnection → Deprotection → Coupling
    """
    # Initialize flags for each transformation
    found_reduction = False
    found_disconnection = False
    found_deprotection = False
    found_coupling = False

    # Initialize depths for each transformation
    reduction_depth = -1
    disconnection_depth = -1
    deprotection_depth = -1
    coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_reduction, found_disconnection, found_deprotection, found_coupling
        nonlocal reduction_depth, disconnection_depth, deprotection_depth, coupling_depth

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)

                # Check for reduction (C=O to CH2-OH)
                carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")
                alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        product_mol
                        and reactant_mol
                        and product_mol.HasSubstructMatch(carbonyl_pattern)
                        and reactant_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        found_reduction = True
                        reduction_depth = depth
                        print(f"Found reduction at depth {depth}")

                # Check for alkene disconnection
                alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
                alkyne_pattern = Chem.MolFromSmarts("[#6]#[#6]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        product_mol
                        and reactant_mol
                        and product_mol.HasSubstructMatch(alkene_pattern)
                        and reactant_mol.HasSubstructMatch(alkyne_pattern)
                    ):
                        found_disconnection = True
                        disconnection_depth = depth
                        print(f"Found disconnection at depth {depth}")

                # Check for ester hydrolysis (deprotection)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=O)[OH]")
                methyl_ester_pattern = Chem.MolFromSmarts("[#6](=O)[O][CH3]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        product_mol
                        and reactant_mol
                        and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                        and reactant_mol.HasSubstructMatch(methyl_ester_pattern)
                    ):
                        found_deprotection = True
                        deprotection_depth = depth
                        print(f"Found deprotection at depth {depth}")

                # Check for sulfonamide coupling
                if len(reactants) > 1:
                    sulfonamide_pattern = Chem.MolFromSmarts("[NH2][S](=O)(=O)[#6]")

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                            found_coupling = True
                            coupling_depth = depth
                            print(f"Found coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if all transformations were found and in the correct sequence
    correct_sequence = (
        found_reduction
        and found_disconnection
        and found_deprotection
        and found_coupling
        and reduction_depth > disconnection_depth > deprotection_depth > coupling_depth
    )

    return correct_sequence
