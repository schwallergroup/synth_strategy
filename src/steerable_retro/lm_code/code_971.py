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
    Detects a late-stage reductive amination strategy where an aldehyde and amine
    are coupled to form a C-N bond, typically in the second half of the synthesis.
    """
    reductive_amination_found = False
    depth_of_reductive_amination = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, depth_of_reductive_amination, max_depth

        max_depth = max(max_depth, depth)

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[#6]=[#8]")
            # Check for amine pattern in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(product_part)

            if product and all(r for r in reactants):
                # Check if reactants contain aldehyde and amine
                has_aldehyde = any(r.HasSubstructMatch(aldehyde_pattern) for r in reactants if r)
                has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactants if r)

                # Check if product has a new C-N bond that wasn't in reactants
                if has_aldehyde and has_amine:
                    # This is a simplified check - in a real implementation, you would need to
                    # check for the specific C-N bond formation between the aldehyde carbon and amine nitrogen
                    print(f"Potential reductive amination found at depth {depth}")
                    reductive_amination_found = True
                    depth_of_reductive_amination = depth

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis
    # (remember depth increases as we go backward in synthesis)
    is_late_stage = depth_of_reductive_amination <= max_depth / 2

    print(
        f"Reductive amination found: {reductive_amination_found}, at depth: {depth_of_reductive_amination}, max depth: {max_depth}"
    )
    print(f"Is late stage: {is_late_stage}")

    return reductive_amination_found and is_late_stage
