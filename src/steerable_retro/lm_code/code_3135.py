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
    This function detects if the synthetic route involves heterocycle formation
    after an amide coupling step.
    """
    heterocycle_formation = False
    amide_coupling = False
    heterocycle_depth = -1
    amide_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation, amide_coupling, heterocycle_depth, amide_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for oxazole formation
                    oxazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#8]1")
                    if product_mol.HasSubstructMatch(oxazole_pattern):
                        # Check if reactants don't have the oxazole pattern
                        oxazole_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(oxazole_pattern):
                                oxazole_in_reactants = True

                        if not oxazole_in_reactants:
                            print(f"Heterocycle formation detected at depth {depth}")
                            heterocycle_formation = True
                            heterocycle_depth = depth

                    # Check for amide coupling
                    amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        # Check if reactants contain carboxylic acid and amine
                        has_acid = False
                        has_amine = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8]")
                                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                                if reactant_mol.HasSubstructMatch(acid_pattern):
                                    has_acid = True
                                if reactant_mol.HasSubstructMatch(amine_pattern):
                                    has_amine = True

                        if has_acid and has_amine:
                            print(f"Amide coupling detected at depth {depth}")
                            amide_coupling = True
                            amide_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if heterocycle formation occurs after amide coupling
    result = heterocycle_formation and amide_coupling and heterocycle_depth < amide_depth

    print(f"Heterocycle formation after amide coupling: {result}")
    return result
