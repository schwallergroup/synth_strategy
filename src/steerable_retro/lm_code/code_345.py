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
    Detects if the synthetic route involves a convergent approach where two complex fragments
    are joined together, rather than a purely linear synthesis.
    """
    found_convergent_step = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If a reaction has multiple reactants, it might be convergent
            if len(reactants) >= 2:
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Filter out None molecules
                    valid_mols = [mol for mol in reactant_mols if mol is not None]

                    # Check complexity of each fragment (using number of atoms as a simple metric)
                    complex_fragments = 0
                    for mol in valid_mols:
                        if (
                            mol.GetNumAtoms() >= 6
                        ):  # Consider fragments with at least 6 atoms as complex
                            complex_fragments += 1

                    if complex_fragments >= 2:
                        print(f"Found convergent assembly of complex fragments at depth {depth}")
                        found_convergent_step = True
                except:
                    print("Error processing reaction SMILES")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_convergent_step
