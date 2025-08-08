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
    This function detects alkylation of a heterocyclic core structure.
    Specifically, it looks for addition of an alkyl chain to a nitrogen-containing aromatic ring.
    """
    alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if an alkyl halide is present in reactants
                alkyl_halide_pattern = re.compile(r".*[CI][CH2][CH2].*")

                # Check if product has a heterocycle
                heterocycle_pattern = Chem.MolFromSmarts("[n;r5,r6]")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(heterocycle_pattern):
                        for reactant in reactants:
                            if alkyl_halide_pattern.match(reactant):
                                print(f"Heterocycle alkylation detected: {rsmi}")
                                alkylation_detected = True
                                break
                except:
                    print("Error processing SMILES in alkylation detection")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return alkylation_detected
