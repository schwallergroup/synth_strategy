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
    Detects if the synthesis route involves an amide bond formation step.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, the product has the amide bond that gets broken
                amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
                carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
                amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                        # Check if reactants contain carboxylic acid and amine
                        has_acid = False
                        has_amine = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if not reactant_mol:
                                continue

                            if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                has_acid = True
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                        if has_acid and has_amine:
                            print(f"Found amide bond formation at depth {depth}")
                            found_amide_formation = True
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_amide_formation
