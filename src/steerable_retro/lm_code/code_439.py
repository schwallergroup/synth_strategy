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
    This function detects if the synthesis involves formation of a quinazolinedione ring system.
    """
    has_quinazolinedione_formation = False

    def dfs_traverse(node):
        nonlocal has_quinazolinedione_formation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains quinazolinedione pattern
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                quinazolinedione_pattern = Chem.MolFromSmarts(
                    "[#6]1[#7][#6](=[O])[#7][#6](=[O])[#6]1"
                )

                # Check if any reactant contains the pattern
                reactant_has_pattern = False
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(quinazolinedione_pattern):
                        reactant_has_pattern = True
                        break

                # If product has pattern but reactants don't, it's a formation
                if (
                    product_mol.HasSubstructMatch(quinazolinedione_pattern)
                    and not reactant_has_pattern
                ):
                    has_quinazolinedione_formation = True
                    print(f"Found quinazolinedione formation: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_quinazolinedione_formation
