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
    This function detects if the synthesis route involves a cross-coupling reaction
    in the early stages (higher depth values).
    """
    cross_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal cross_coupling_detected

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth >= 2)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid in reactants (indicator of Suzuki coupling)
                boronic_acid_pattern = re.compile(r"B\(O\)|OB\(O\)")
                has_boronic_acid = any(
                    boronic_acid_pattern.search(reactant) for reactant in reactants
                )

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I]")
                has_aryl_halide = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                        break

                if has_boronic_acid and has_aryl_halide:
                    cross_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Cross-coupling in early stage detected: {cross_coupling_detected}")
    return cross_coupling_detected
