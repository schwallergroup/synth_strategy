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
    This function detects SNAr reaction where an amine displaces a halide.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for aryl halide pattern
                aryl_halide_pattern = Chem.MolFromSmarts("c-[Br,Cl,I,F]")
                # Check for amine pattern
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                # Check for aryl amine pattern in product
                aryl_amine_pattern = Chem.MolFromSmarts("c-[NH]")

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
                    print("Detected SNAr reaction with amine")
                    snar_found = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return snar_found
