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
    This function detects if the synthetic route involves the strategic use of nitrile groups.
    """
    nitrile_present = False
    nitrile_transformed = False

    def dfs_traverse(node):
        nonlocal nitrile_present, nitrile_transformed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            for reactant in reactants_smiles:
                if reactant.strip():
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_present = True
                        print(f"Nitrile group found in reactant: {reactant}")

            # Check if nitrile is transformed
            if product_smiles and nitrile_present:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # Check if product has different nitrile pattern or count than reactants
                    product_nitrile_count = len(product_mol.GetSubstructMatches(nitrile_pattern))

                    reactants_nitrile_count = 0
                    for reactant in reactants_smiles:
                        if reactant.strip():
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                reactants_nitrile_count += len(
                                    reactant_mol.GetSubstructMatches(nitrile_pattern)
                                )

                    if product_nitrile_count != reactants_nitrile_count:
                        nitrile_transformed = True
                        print(
                            f"Nitrile transformation detected: Reactants have {reactants_nitrile_count} nitriles, product has {product_nitrile_count}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both conditions are met
    result = nitrile_present and nitrile_transformed
    print(f"Nitrile utilization strategy detected: {result}")
    return result
