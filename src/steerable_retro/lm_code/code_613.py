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
    This function detects a strategy involving late-stage amide formation,
    specifically looking for amide bond formation in the final steps of the synthesis.
    """
    has_late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]-[#6]")
                acyl_pattern = Chem.MolFromSmarts("[#6](=[O])-[Cl,Br,I,O]")
                amide_pattern = Chem.MolFromSmarts("[#7]-[#6](=[O])")

                # Check reactants for amine and acyl groups
                has_amine = False
                has_acyl = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if mol.HasSubstructMatch(acyl_pattern):
                                has_acyl = True
                    except:
                        continue

                # Check product for amide group
                product_has_amide = False
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(amide_pattern):
                        product_has_amide = True
                except:
                    pass

                # If amide is formed and it's in the late stage (depth <= 1)
                if has_amine and has_acyl and product_has_amide and depth <= 1:
                    has_late_amide_formation = True
                    print(f"Detected late-stage amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_amide_formation
