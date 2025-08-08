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
    Detects convergent synthesis with late-stage N-alkylation joining two complex fragments,
    where one fragment contains a piperazine and the other contains a chloroalkyl group.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track if we found N-alkylation
    n_alkylation_found = False
    # Track if we found piperazine fragment
    piperazine_fragment_found = False
    # Track if we found chloroalkyl fragment
    chloroalkyl_fragment_found = False

    def dfs_traverse(node):
        nonlocal found_pattern, n_alkylation_found, piperazine_fragment_found, chloroalkyl_fragment_found

        if node["type"] == "reaction":
            # Check if this is a reaction node with metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-alkylation reaction
                # Pattern: secondary amine + chloroalkyl -> tertiary amine
                if len(reactants) >= 2:
                    # Check for piperazine in reactants
                    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
                    # Check for chloroalkyl in reactants
                    chloroalkyl_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#17]")

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                if mol.HasSubstructMatch(piperazine_pattern):
                                    piperazine_fragment_found = True
                                if mol.HasSubstructMatch(chloroalkyl_pattern):
                                    chloroalkyl_fragment_found = True
                        except:
                            continue

                    # Check if product has a tertiary amine formed from joining these fragments
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            # Check if product has both piperazine and former chloroalkyl (now connected)
                            if (
                                prod_mol.HasSubstructMatch(piperazine_pattern)
                                and piperazine_fragment_found
                                and chloroalkyl_fragment_found
                            ):
                                n_alkylation_found = True
                                print(
                                    "Found N-alkylation joining piperazine and chloroalkyl fragments"
                                )
                    except:
                        pass

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found the complete pattern
    found_pattern = n_alkylation_found and piperazine_fragment_found and chloroalkyl_fragment_found

    if found_pattern:
        print("Detected convergent synthesis with late-stage N-alkylation")

    return found_pattern
