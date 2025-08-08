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
    This function detects a synthetic strategy where a nitro group is carried through
    multiple steps and reduced to an amine in the final or near-final step.
    """
    nitro_reduction_depth = None
    max_depth = -1
    nitro_present = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, max_depth, nitro_present

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            reactants_with_nitro = 0
            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[#7+](=[#8])[#8-]")):
                        reactants_with_nitro += 1
                        nitro_present = True
                except:
                    continue

            # Check for amine in product but not in reactants (nitro reduction)
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                if p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                    # Check if reactants had nitro but not amine
                    if reactants_with_nitro > 0:
                        reactants_with_amine = 0
                        for r_smi in reactants_smiles:
                            r_mol = Chem.MolFromSmiles(r_smi)
                            if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                                reactants_with_amine += 1

                        if reactants_with_amine == 0:
                            nitro_reduction_depth = depth
                            print(f"Nitro reduction detected at depth {depth}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro reduction occurred in the final or near-final step
    if nitro_reduction_depth is not None and nitro_present:
        # Consider it late-stage if it's in the first third of the synthesis depth
        if nitro_reduction_depth <= max_depth / 3:
            print(
                f"Late-stage nitro reduction strategy detected (depth {nitro_reduction_depth} of {max_depth})"
            )
            return True

    return False
