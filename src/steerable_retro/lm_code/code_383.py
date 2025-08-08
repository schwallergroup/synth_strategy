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
    This function detects a synthetic strategy involving heterocycle elaboration
    through Suzuki coupling (aryl-heteroaryl C-C bond formation).
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Suzuki coupling pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for boronic acid/ester in reactants
            boronic_pattern = Chem.MolFromSmarts("[#6]B([O])[O]")
            boronic_ester_pattern = Chem.MolFromSmarts("[#6]B1OC(C)(C)OC1(C)C")

            # Look for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

            # Look for heterocycle pattern (imidazole)
            imidazole_pattern = Chem.MolFromSmarts("c1ncnc1")

            has_boronic = any(
                mol
                and (
                    mol.HasSubstructMatch(boronic_pattern)
                    or mol.HasSubstructMatch(boronic_ester_pattern)
                )
                for mol in reactant_mols
            )

            has_aryl_halide = any(
                mol and mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols
            )

            has_imidazole = any(
                mol and mol.HasSubstructMatch(imidazole_pattern) for mol in reactant_mols
            )

            if has_boronic and has_aryl_halide and has_imidazole:
                print("Detected Suzuki coupling for heterocycle elaboration")
                has_suzuki_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_suzuki_coupling
