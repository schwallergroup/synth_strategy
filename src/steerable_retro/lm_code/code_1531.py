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
    Detects if the synthesis route involves a sequence of functional group
    interconversions: ester → acid → amide.
    """
    # Track functional groups at each depth
    functional_groups = {}

    def dfs_traverse(node, depth=0):
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check for functional groups
                ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")
                acid_pattern = Chem.MolFromSmarts("[C](=O)[O;H,-]")
                amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")

                has_ester = product_mol.HasSubstructMatch(ester_pattern) if ester_pattern else False
                has_acid = product_mol.HasSubstructMatch(acid_pattern) if acid_pattern else False
                has_amide = product_mol.HasSubstructMatch(amide_pattern) if amide_pattern else False

                functional_groups[depth] = {
                    "ester": has_ester,
                    "acid": has_acid,
                    "amide": has_amide,
                }

                if has_ester:
                    print(f"Ester found at depth {depth}")
                if has_acid:
                    print(f"Carboxylic acid found at depth {depth}")
                if has_amide:
                    print(f"Amide found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check for the sequence: ester → acid → amide
    # Note: In retrosynthetic direction (increasing depth), this would be amide → acid → ester
    depths = sorted(functional_groups.keys())

    sequence_found = False
    for i in range(len(depths) - 2):
        d1, d2, d3 = depths[i], depths[i + 1], depths[i + 2]

        if (
            functional_groups[d1]["amide"]
            and functional_groups[d2]["acid"]
            and functional_groups[d3]["ester"]
        ):
            sequence_found = True
            print(f"Found ester→acid→amide sequence at depths {d3}→{d2}→{d1}")
            break

    return sequence_found
