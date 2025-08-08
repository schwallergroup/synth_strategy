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
    This function detects if the synthetic route employs late-stage aromatic functionalization.
    It looks for reactions in the first half of the synthesis that modify aromatic rings.
    """
    max_depth = 0
    aromatic_func_depths = []

    # First pass to determine max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find aromatic functionalizations
    def find_aromatic_func(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic patterns
            aromatic_pattern = Chem.MolFromSmarts("[c]")

            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(aromatic_pattern):
                    # Check if this reaction modifies an aromatic ring
                    for reactant in reactants:
                        try:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and react_mol.HasSubstructMatch(aromatic_pattern):
                                # This is a reaction involving aromatic rings
                                aromatic_func_depths.append(depth)
                                print(f"Found aromatic functionalization at depth {depth}: {rsmi}")
                                break
                        except:
                            continue
            except:
                pass

        for child in node.get("children", []):
            find_aromatic_func(child, depth + 1)

    find_max_depth(route)
    find_aromatic_func(route)

    # Check if there are aromatic functionalizations in the first half of the synthesis
    if aromatic_func_depths:
        half_depth = max_depth / 2
        late_stage_funcs = [d for d in aromatic_func_depths if d < half_depth]
        if late_stage_funcs:
            print(f"Found late-stage aromatic functionalizations at depths: {late_stage_funcs}")
            return True

    return False
