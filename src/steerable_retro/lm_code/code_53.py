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
    This function detects convergent synthesis with multiple fragment couplings.
    """
    coupling_reactions = 0

    def dfs_traverse(node):
        nonlocal coupling_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, it might be a coupling reaction
            if len(reactants) >= 2:
                # Check for common coupling patterns
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")
                boronic_acid_pattern = Chem.MolFromSmarts("[c]B(O)O")
                carbonyl_pattern = Chem.MolFromSmarts("[C]=O")

                has_aryl_halide = False
                has_boronic_acid = False
                has_carbonyl = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if mol.HasSubstructMatch(boronic_acid_pattern):
                                has_boronic_acid = True
                            if mol.HasSubstructMatch(carbonyl_pattern):
                                has_carbonyl = True
                    except:
                        continue

                # Identify coupling reactions
                if (has_aryl_halide and has_boronic_acid) or (has_carbonyl and len(reactants) >= 2):
                    coupling_reactions += 1
                    print(f"Coupling reaction detected, total: {coupling_reactions}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If there are multiple coupling reactions, it's a convergent synthesis
    if coupling_reactions >= 2:
        print("Convergent synthesis strategy with multiple fragment couplings detected")
        return True
    return False
