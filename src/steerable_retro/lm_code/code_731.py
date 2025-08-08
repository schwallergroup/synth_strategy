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
    This function detects multiple heteroatom functionalization steps in the synthesis,
    including sulfonamide formation, N-methylation, and amide formation.
    """
    functionalization_count = 0

    def dfs_traverse(node):
        nonlocal functionalization_count

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonamide formation
            sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])")
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[Cl]")

            # Check for N-methylation
            n_methyl_pattern = Chem.MolFromSmarts("[#7]-[#6]")

            # Check for amide formation
            amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])")
            acyl_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])-[Cl]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Check for sulfonamide formation
                if product_mol.HasSubstructMatch(sulfonamide_pattern):
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                            functionalization_count += 1
                            print(f"Detected sulfonamide formation: {rsmi}")
                            break

                # Check for N-methylation
                if product_mol.HasSubstructMatch(n_methyl_pattern):
                    methyl_iodide_found = False
                    amine_found = False
                    for reactant in reactants:
                        if "CI" in reactant or "IC" in reactant:
                            methyl_iodide_found = True
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[#7;H1,H2]")):
                            amine_found = True

                    if methyl_iodide_found and amine_found:
                        functionalization_count += 1
                        print(f"Detected N-methylation: {rsmi}")

                # Check for amide formation
                if product_mol.HasSubstructMatch(amide_pattern):
                    acyl_chloride_found = False
                    amine_found = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(acyl_chloride_pattern):
                                acyl_chloride_found = True
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#7;H1,H2]")):
                                amine_found = True

                    if acyl_chloride_found and amine_found:
                        functionalization_count += 1
                        print(f"Detected amide formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return functionalization_count >= 2
