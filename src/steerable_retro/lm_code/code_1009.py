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
    This function detects if the synthetic route involves nucleophilic aromatic substitution
    of chloro to amino on a pyrimidine ring.
    """
    nas_detected = False

    def dfs_traverse(node):
        nonlocal nas_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain chloropyrimidine
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                reactant_with_chloro = None
                for r in reactants:
                    if r is not None and r.HasSubstructMatch(Chem.MolFromSmarts("[#17]-c1ncccn1")):
                        reactant_with_chloro = r
                        break

                # Check if product contains aminopyrimidine
                product = Chem.MolFromSmiles(product_smiles)
                if reactant_with_chloro is not None and product is not None:
                    if product.HasSubstructMatch(Chem.MolFromSmarts("[#7H2]-c1ncccn1")):
                        nas_detected = True
                        print("Nucleophilic aromatic substitution detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nas_detected
