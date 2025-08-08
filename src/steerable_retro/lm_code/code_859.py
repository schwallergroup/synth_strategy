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
    This function detects if the synthetic route involves aromatic C-N coupling reactions.
    """
    aromatic_cn_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_cn_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic C-N coupling pattern
                # Look for aryl halide and aniline patterns
                aryl_halide_pattern = Chem.MolFromSmarts("c[Cl,Br,I,F]")
                aniline_pattern = Chem.MolFromSmarts("c[NH2]")

                product_mol = Chem.MolFromSmiles(product)
                diarylamine_pattern = Chem.MolFromSmarts("c[NH]c")

                if product_mol and product_mol.HasSubstructMatch(diarylamine_pattern):
                    has_aryl_halide = False
                    has_aniline = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                has_aryl_halide = True
                            if reactant_mol.HasSubstructMatch(aniline_pattern):
                                has_aniline = True

                    if has_aryl_halide and has_aniline:
                        aromatic_cn_coupling_found = True
                        print(f"Aromatic C-N coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return aromatic_cn_coupling_found
