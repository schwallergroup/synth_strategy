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
    Detects if the synthesis route employs a ketone formation strategy via acylation,
    where a C-C bond is formed creating a ketone between aromatic rings.
    """
    ketone_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ketone_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O;H1]")
                # Check for ketone between aromatics in product
                ketone_pattern = Chem.MolFromSmarts("[#6;a][C](=[O])[#6;a]")

                carboxylic_acid_present = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                        carboxylic_acid_present = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    carboxylic_acid_present
                    and product_mol
                    and product_mol.HasSubstructMatch(ketone_pattern)
                ):
                    print(f"Ketone formation via acylation detected at depth {depth}")
                    ketone_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ketone_formation_detected
