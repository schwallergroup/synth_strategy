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
    This function detects if there are sequential SNAr reactions (aryl chloride displacements).
    """
    snar_reactions = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl chloride in reactants
                aryl_cl_pattern = Chem.MolFromSmarts("c-Cl")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                has_aryl_cl = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_cl_pattern):
                            has_aryl_cl = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # Check if product has new C-N bond where chlorine was
                if has_aryl_cl and has_amine:
                    snar_reactions.append(node.get("metadata", {}).get("ID", "unknown"))

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if there are consecutive SNAr reactions
    print(f"Found {len(snar_reactions)} potential SNAr reactions")
    return len(snar_reactions) >= 2
