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
    This function detects an aryl-amine coupling (Buchwald-Hartwig type) in the synthesis route.
    """
    aryl_amine_coupling_found = False

    def dfs_traverse(node):
        nonlocal aryl_amine_coupling_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants include aryl halide and amine, and product has aryl-amine bond
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I]")
                amine_pattern = Chem.MolFromSmarts("[N;H2]")
                aryl_amine_pattern = Chem.MolFromSmarts("[c]-[N;H]")

                reactant_has_aryl_halide = any(
                    mol and mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols
                )
                reactant_has_amine = any(
                    mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                )
                product_has_aryl_amine = product_mol and product_mol.HasSubstructMatch(
                    aryl_amine_pattern
                )

                if reactant_has_aryl_halide and reactant_has_amine and product_has_aryl_amine:
                    print("Aryl-amine coupling detected")
                    aryl_amine_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aryl_amine_coupling_found
