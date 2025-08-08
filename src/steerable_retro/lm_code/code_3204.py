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
    Detects a strategy involving late-stage disconnection of an ether bond between aromatic fragments.
    Looks for an ether cleavage reaction at a shallow depth (late in the synthesis).
    """
    found_late_stage_ether_disconnection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_ether_disconnection

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ether bond between aromatic systems
                ether_pattern = Chem.MolFromSmarts("[c]-[O]-[CH2]-[c]")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check if product has the ether bond but at least one reactant doesn't
                    if (
                        product_mol
                        and product_mol.HasSubstructMatch(ether_pattern)
                        and any(
                            mol and not mol.HasSubstructMatch(ether_pattern)
                            for mol in reactant_mols
                        )
                    ):
                        found_late_stage_ether_disconnection = True
                        print(f"Found late-stage ether disconnection at depth {depth}: {rsmi}")
                except:
                    pass

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_stage_ether_disconnection
