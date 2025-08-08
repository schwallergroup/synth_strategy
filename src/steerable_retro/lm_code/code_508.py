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
    Detects a synthetic strategy involving aromatic nucleophilic substitution
    with a piperazine nucleophile.
    """
    aromatic_snar_detected = False

    def dfs_traverse(node):
        nonlocal aromatic_snar_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloropyridine in reactants
                chloropyridine_present = False
                piperazine_present = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        chloropyridine_pattern = Chem.MolFromSmarts(
                            "[Cl][#6]1[#6][#6][#6][#7][#6]1"
                        )
                        piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

                        if reactant_mol.HasSubstructMatch(chloropyridine_pattern):
                            chloropyridine_present = True
                        if reactant_mol.HasSubstructMatch(piperazine_pattern):
                            piperazine_present = True

                # Check for piperazine-linked pyridine in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    linked_pattern = Chem.MolFromSmarts(
                        "[#6]1[#6][#6][#6][#7][#6]1[#7]2[#6][#6][#7][#6][#6]2"
                    )
                    if (
                        product_mol.HasSubstructMatch(linked_pattern)
                        and chloropyridine_present
                        and piperazine_present
                    ):
                        aromatic_snar_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Aromatic nucleophilic substitution detected: {aromatic_snar_detected}")

    return aromatic_snar_detected
