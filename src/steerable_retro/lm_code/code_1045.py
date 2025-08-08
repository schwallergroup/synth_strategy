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
    Detects if the synthetic route involves a late-stage Suzuki coupling (in the final
    or penultimate step) to introduce complex aryl fragments.
    """
    borylation_found = False
    suzuki_found = False
    borylation_depth = -1
    suzuki_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal borylation_found, suzuki_found, borylation_depth, suzuki_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for borylation (aryl halide to boronic acid/ester)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]-[B]([O])[O]")
                ):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]-[Br,Cl,I]")
                        ):
                            print(f"Found borylation at depth {depth}")
                            borylation_found = True
                            borylation_depth = depth

                # Check for Suzuki coupling (boronic acid/ester + aryl halide -> biaryl)
                has_boronic = False
                has_aryl_halide = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[B]([O])[O]")):
                            has_boronic = True
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[Br,Cl,I]")):
                            has_aryl_halide = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_boronic
                    and has_aryl_halide
                    and product_mol
                    and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[#6]"))
                ):
                    print(f"Found Suzuki coupling at depth {depth}")
                    suzuki_found = True
                    suzuki_depth = depth

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if Suzuki is in the final (depth 0) or penultimate step (depth 1)
    is_late_stage = suzuki_depth <= 1

    # Check if borylation precedes Suzuki
    correct_sequence = borylation_depth > suzuki_depth

    return suzuki_found and borylation_found and is_late_stage and correct_sequence
