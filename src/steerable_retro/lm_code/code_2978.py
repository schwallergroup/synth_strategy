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
    Detects a synthetic strategy involving Boc protection of an amine.
    """
    has_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check if this is a Boc protection reaction
            boc_protection_reactions = [
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of secondary amine",
                "Boc amine protection of primary amine",
            ]

            for reaction_name in boc_protection_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    print(f"Detected {reaction_name} reaction at depth {depth}: {rsmi}")
                    has_boc_protection = True
                    break

            # If we didn't find a specific Boc protection reaction, check for Boc group addition
            if not has_boc_protection:
                try:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if any reactant has an amine and product has a Boc group
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and checker.check_fg("Boc", product_smiles):
                        # Check if any reactant has a primary or secondary amine
                        for reactant in reactants_smiles:
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                print(
                                    f"Detected Boc protection (amine to Boc) at depth {depth}: {rsmi}"
                                )
                                has_boc_protection = True
                                break
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_boc_protection
