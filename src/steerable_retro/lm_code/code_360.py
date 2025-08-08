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
    Detects if the synthesis route involves coupling of a fragment with an amidine derivative.
    """
    amidine_coupling = False

    def dfs_traverse(node):
        nonlocal amidine_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Need at least 2 reactants for fragment coupling
                if len(reactants) >= 2:
                    # Check for amidine in reactants
                    amidine_pattern = Chem.MolFromSmarts("[#6](=[#7])-[#7]")
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(amidine_pattern):
                                # Check if another reactant has a sulfonamide
                                sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])")
                                for other_reactant in reactants:
                                    if other_reactant != reactant:
                                        try:
                                            other_mol = Chem.MolFromSmiles(other_reactant)
                                            if other_mol and other_mol.HasSubstructMatch(
                                                sulfonamide_pattern
                                            ):
                                                amidine_coupling = True
                                        except:
                                            continue
                        except:
                            continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Fragment coupling with amidine: {amidine_coupling}")
    return amidine_coupling
