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
    Detects a strategy involving multiple sequential dehalogenation steps
    (removal of halogens like Br, I).
    """
    dehalogenation_count = 0

    def dfs_traverse(node):
        nonlocal dehalogenation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Find the main reactant (usually the first one)
                main_reactant = None
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.GetNumHeavyAtoms() > 5:  # Assuming main reactant is larger
                        main_reactant = r_mol
                        break

                if product_mol and main_reactant:
                    # Check for halogen removal
                    br_pattern = Chem.MolFromSmarts("[#6]-[Br]")
                    i_pattern = Chem.MolFromSmarts("[#6]-[I]")

                    # Check if reactant has halogen that product doesn't have
                    if (
                        main_reactant.HasSubstructMatch(br_pattern)
                        and not product_mol.HasSubstructMatch(br_pattern)
                    ) or (
                        main_reactant.HasSubstructMatch(i_pattern)
                        and not product_mol.HasSubstructMatch(i_pattern)
                    ):
                        dehalogenation_count += 1
                        print(f"Detected dehalogenation step, total count: {dehalogenation_count}")

            except Exception as e:
                print(f"Error in SMILES processing: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return dehalogenation_count >= 2  # Return True if at least 2 dehalogenation steps
