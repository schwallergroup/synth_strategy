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
    Detects if the synthesis route employs a late-stage sulfonylation strategy,
    where a phenol is converted to a sulfonate ester in one of the final steps.
    """
    sulfonylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonylation_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Focus on late-stage reactions (depth 0 or 1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[#6][O;H1]")
                # Check for sulfonate in product
                sulfonate_pattern = Chem.MolFromSmarts("[#6][O][S](=[O])(=[O])[#6]")

                reactant_has_phenol = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phenol_pattern):
                        reactant_has_phenol = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_phenol
                    and product_mol
                    and product_mol.HasSubstructMatch(sulfonate_pattern)
                ):
                    print(f"Late-stage sulfonylation detected at depth {depth}")
                    sulfonylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return sulfonylation_detected
