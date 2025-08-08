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
    Detects a strategy involving ketone reduction to alcohol.
    """
    has_ketone_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ketone_reduction

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                # Check if this is a ketone reduction reaction
                if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                    print(f"Detected ketone reduction reaction at depth {depth}")
                    has_ketone_reduction = True
                else:
                    # Alternative check using functional groups
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                    product_mol = Chem.MolFromSmiles(products_part)

                    if all(reactant_mols) and product_mol:
                        # Check for ketone in reactants
                        has_ketone_in_reactants = any(
                            checker.check_fg("Ketone", Chem.MolToSmiles(mol))
                            for mol in reactant_mols
                        )

                        # Check for alcohol in product
                        has_alcohol_in_product = checker.check_fg(
                            "Secondary alcohol", products_part
                        )

                        if has_ketone_in_reactants and has_alcohol_in_product:
                            # Additional check to ensure it's not a false positive
                            # This is a reduction reaction converting ketone to alcohol
                            print(
                                f"Detected ketone reduction via functional group analysis at depth {depth}"
                            )
                            has_ketone_reduction = True
            except Exception as e:
                print(f"Error analyzing ketone reduction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Ketone reduction strategy detected: {has_ketone_reduction}")
    return has_ketone_reduction
