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
    This function detects if the synthesis includes coupling reactions involving heteroaromatic rings,
    particularly pyridine-containing compounds.
    """
    heteroaromatic_coupling_detected = False

    def dfs_traverse(node):
        nonlocal heteroaromatic_coupling_detected

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for heteroaromatic rings in reactants
            pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")  # Simple pyridine pattern

            # Check for biaryl formation
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_heteroaromatic = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(pyridine_pattern):
                        has_heteroaromatic = True
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[#53,#35]")):
                        has_aryl_halide = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_heteroaromatic
                    and has_aryl_halide
                    and product_mol
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    print("Detected heteroaromatic coupling reaction")
                    heteroaromatic_coupling_detected = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return heteroaromatic_coupling_detected
