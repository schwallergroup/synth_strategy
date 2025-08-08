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
    This function detects late-stage reductive amination.
    Looks for C-N bond formation between an aldehyde and a secondary amine
    in the final steps of the synthesis.
    """
    reductive_amination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")

                # Check for amine in reactants (focusing on piperazine)
                amine_pattern = Chem.MolFromSmarts("N1CCN(C)CC1")

                # Check for C-N bond in product
                cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                has_aldehyde = False
                has_amine = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(aldehyde_pattern):
                            has_aldehyde = True
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    has_cn_bond = product_mol and product_mol.HasSubstructMatch(cn_bond_pattern)

                    if has_aldehyde and has_amine and has_cn_bond:
                        print("Detected late-stage reductive amination")
                        reductive_amination_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return reductive_amination_detected
