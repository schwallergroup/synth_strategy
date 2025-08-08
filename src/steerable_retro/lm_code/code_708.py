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
    This function detects the incorporation of fluorinated fragments (CF3, CF2)
    in the late stages of synthesis (lower depth in the tree).
    """
    fluorinated_fragments_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal fluorinated_fragments_late_stage

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Consider only late-stage reactions (depth 0, 1, 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant contains fluorine
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    for reactant in reactants:
                        if reactant:
                            # Check for CF3 or CF2 groups
                            cf3_pattern = Chem.MolFromSmarts("[#6]([F])([F])[F]")
                            cf2_pattern = Chem.MolFromSmarts("[#6]([F])([F])")

                            if reactant.HasSubstructMatch(
                                cf3_pattern
                            ) or reactant.HasSubstructMatch(cf2_pattern):
                                print(
                                    f"Detected fluorinated fragment in late-stage reaction at depth {depth}"
                                )
                                fluorinated_fragments_late_stage = True
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return fluorinated_fragments_late_stage
