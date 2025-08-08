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
    This function detects if the synthesis includes a sequence of benzylic functionalization
    steps (e.g., CH3 → CH2Br → CH2NH2).
    """
    benzylic_steps = []

    def dfs_traverse(node):
        nonlocal benzylic_steps

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Patterns for benzylic transformations
            benzylic_methyl_pattern = Chem.MolFromSmarts("[c][C]")
            benzylic_bromomethyl_pattern = Chem.MolFromSmarts("[c][C][Br]")
            benzylic_aminomethyl_pattern = Chem.MolFromSmarts("[c][C][N]")

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                # Check for CH3 → CH2Br transformation
                if (
                    any(
                        r and r.HasSubstructMatch(benzylic_methyl_pattern)
                        for r in reactant_mols
                        if r
                    )
                    and product_mol
                    and product_mol.HasSubstructMatch(benzylic_bromomethyl_pattern)
                ):
                    benzylic_steps.append("methyl_to_bromomethyl")
                    print("Detected benzylic methyl to bromomethyl transformation")

                # Check for CH2Br → CH2NH2 transformation
                if (
                    any(
                        r and r.HasSubstructMatch(benzylic_bromomethyl_pattern)
                        for r in reactant_mols
                        if r
                    )
                    and product_mol
                    and product_mol.HasSubstructMatch(benzylic_aminomethyl_pattern)
                ):
                    benzylic_steps.append("bromomethyl_to_aminomethyl")
                    print("Detected benzylic bromomethyl to aminomethyl transformation")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have the complete sequence
    return (
        "methyl_to_bromomethyl" in benzylic_steps and "bromomethyl_to_aminomethyl" in benzylic_steps
    )
