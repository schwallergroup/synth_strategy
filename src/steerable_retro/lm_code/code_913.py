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
    This function detects if there's a late-stage amide coupling introducing a nitrile-containing side chain.
    """
    found_late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amide_coupling

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider reactions at depth 0 or 1 (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#8;H1]")
                nitrile_pattern = Chem.MolFromSmarts("[#6]-C#N")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#7]")

                has_acid = False
                has_nitrile = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(acid_pattern):
                            has_acid = True
                        if mol.HasSubstructMatch(nitrile_pattern):
                            has_nitrile = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # Check if product has amide bond
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(amide_pattern)
                    and product_mol.HasSubstructMatch(nitrile_pattern)
                ):
                    if has_acid and has_amine and has_nitrile:
                        found_late_amide_coupling = True
                        print("Found late-stage amide coupling with nitrile group")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_late_amide_coupling
