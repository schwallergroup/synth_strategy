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
    Detects a strategy where a phenol directs ortho/para functionalization,
    followed by O-methylation and further transformations.
    """
    has_phenol = False
    has_phenol_functionalization = False
    has_o_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_phenol, has_phenol_functionalization, has_o_methylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in early stages
                if depth >= 4 and not has_phenol:
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and react_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[c]-[OH]")
                        ):
                            has_phenol = True
                            print("Detected phenol at depth", depth)
                            break

                # Check for functionalization of phenol-containing compound
                if depth >= 3 and has_phenol and not has_phenol_functionalization:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and (
                        prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[OH]"))
                        and (
                            prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br]"))
                            or prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[N+](=[O])[O-]"))
                        )
                    ):
                        has_phenol_functionalization = True
                        print("Detected phenol-directed functionalization at depth", depth)

                # Check for O-methylation
                if depth >= 2 and not has_o_methylation:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[O][CH3]")):
                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and react_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[c]-[OH]")
                            ):
                                has_o_methylation = True
                                print("Detected O-methylation at depth", depth)
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have phenol and at least one of the other transformations
    result = has_phenol and (has_phenol_functionalization or has_o_methylation)
    print(f"Phenol-directed functionalization strategy detected: {result}")
    return result
