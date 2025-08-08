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
    Detects a synthetic strategy involving Suzuki coupling to install a vinyl group
    on an aromatic ring.
    """
    suzuki_vinyl_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_vinyl_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aryl bromide in reactants
            aryl_bromide_pattern = Chem.MolFromSmarts("c-Br")
            has_aryl_bromide = False

            # Check for vinyl boronic acid or derivative in reactants
            vinyl_boron_pattern = Chem.MolFromSmarts("C=C-[B]")
            has_vinyl_boron = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_bromide_pattern):
                            has_aryl_bromide = True
                        if (
                            mol.HasSubstructMatch(vinyl_boron_pattern)
                            or "B" in reactant
                            and "C=C" in reactant
                        ):
                            has_vinyl_boron = True
                except:
                    continue

            # Check if product has a vinyl group attached to an aromatic ring
            aryl_vinyl_pattern = Chem.MolFromSmarts("c-C=C")
            has_aryl_vinyl_product = False

            try:
                mol = Chem.MolFromSmiles(product_smiles)
                if mol and mol.HasSubstructMatch(aryl_vinyl_pattern):
                    has_aryl_vinyl_product = True
            except:
                pass

            # If all conditions are met, we have a Suzuki coupling for vinyl installation
            if has_aryl_bromide and has_vinyl_boron and has_aryl_vinyl_product:
                suzuki_vinyl_found = True
                print("Found Suzuki coupling: vinyl group installed on aromatic ring")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return suzuki_vinyl_found
