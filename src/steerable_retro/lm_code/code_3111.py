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
    Detects a convergent synthesis strategy that maintains a CF3 group
    throughout the synthesis.
    """
    found_convergent_step = False

    # Check if target molecule has CF3
    if route["type"] == "mol" and not checker.check_fg("Trifluoro group", route["smiles"]):
        print(f"Target molecule does not have CF3 group")
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for convergent synthesis (multiple fragments coming together)
            if len(reactants_smiles) >= 2:
                # Check if this is a significant fragment coupling
                significant_fragments = 0
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 5
                    ):  # Consider fragments with >5 atoms significant
                        significant_fragments += 1

                if significant_fragments >= 2:
                    found_convergent_step = True
                    print(f"Found convergent synthesis step at depth {depth}")

            # Check if CF3 is maintained in this reaction
            reactant_has_cf3 = any(checker.check_fg("Trifluoro group", r) for r in reactants_smiles)
            product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)

            # If reactant had CF3 but product doesn't, CF3 wasn't maintained
            if reactant_has_cf3 and not product_has_cf3:
                print(f"CF3 group not maintained at depth {depth}")
                return False

        # Continue traversing
        for child in node.get("children", []):
            if not dfs_traverse(child, depth + 1):
                return False

        return True

    # Start traversal and check if CF3 is maintained throughout
    cf3_maintained = dfs_traverse(route)

    return found_convergent_step and cf3_maintained
