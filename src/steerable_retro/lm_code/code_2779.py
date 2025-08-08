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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects if the synthesis follows a linear (non-convergent) strategy.

    A linear synthesis strategy involves sequential transformations of a single main building block,
    while a convergent strategy combines multiple complex fragments at some point.
    """
    is_linear = True

    # List of functional groups commonly found in reagents rather than building blocks
    reagent_fg_types = [
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Magnesium halide",
        "Zinc halide",
        "Tin",
        "Silyl protective group",
        "TMS ether protective group",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, check if they are significant building blocks
            if len(reactants) > 1:
                significant_reactants = []

                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol is None:
                        continue

                    # Skip if it's likely a reagent based on functional groups
                    is_reagent = False
                    for fg in reagent_fg_types:
                        if checker.check_fg(fg, r):
                            is_reagent = True
                            break

                    # Consider as significant if it has enough heavy atoms and isn't a typical reagent
                    if not is_reagent and mol.GetNumHeavyAtoms() > 5:
                        significant_reactants.append(r)

                # If we have multiple significant reactants at a late stage (low depth),
                # this indicates a convergent synthesis
                if len(significant_reactants) > 1 and depth < 3:
                    is_linear = False
                    print(
                        f"Convergent step detected at depth {depth} with {len(significant_reactants)} significant reactants"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
