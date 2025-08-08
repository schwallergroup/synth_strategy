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
    This function detects a strategy involving multiple amide bond formations.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation patterns
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Look for new amide bond in product
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    if product_mol.HasSubstructMatch(amide_pattern):
                        # Check if reactants include acid/acid derivative and amine
                        has_acid = False
                        has_amine = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Check for acid/acid chloride/ester
                                if reactant_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("[C](=[O])[O,Cl]")
                                ):
                                    has_acid = True
                                # Check for amine
                                if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[N]")):
                                    has_amine = True

                        if has_acid and has_amine:
                            amide_formation_count += 1
                            print(f"Detected amide formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if there are multiple amide formations
    if amide_formation_count >= 2:
        print(f"Detected multiple amide formation strategy with {amide_formation_count} instances")
        return True
    return False
