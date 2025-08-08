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
    Detects a strategy involving the incorporation of an aromatic ring
    containing both cyano and chloro substituents.
    """
    cyano_chloro_aromatic_found = False

    def dfs_traverse(node):
        nonlocal cyano_chloro_aromatic_found

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Check for aromatic ring with both cyano and chloro substituents
            cyano_chloro_aromatic_pattern = Chem.MolFromSmarts(
                "[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]#[#7].[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#17]"
            )

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]

            for r in reactants:
                if (
                    r
                    and r.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6]#[#7]")
                    )
                    and r.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#17]")
                    )
                ):
                    print("Cyano and chloro substituted aromatic ring detected")
                    cyano_chloro_aromatic_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Cyano and chloro substituted aromatic ring found: {cyano_chloro_aromatic_found}")

    return cyano_chloro_aromatic_found
