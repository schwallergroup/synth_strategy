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
    This function detects a synthetic strategy involving late-stage deprotection
    of a Boc-protected amine.
    """
    boc_deprotection_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_at_late_stage

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}, RSMI: {rsmi}")

            # Check if this is any type of Boc deprotection reaction
            boc_deprotection_reactions = [
                "Boc amine deprotection",
                "Boc amine deprotection of guanidine",
                "Boc amine deprotection to NH-NH2",
                "Tert-butyl deprotection of amine",
            ]

            is_boc_deprotection = any(
                checker.check_reaction(rxn, rsmi) for rxn in boc_deprotection_reactions
            )
            print(f"Is Boc deprotection: {is_boc_deprotection}")

            # If not detected by reaction type, try to detect by functional group changes
            if not is_boc_deprotection:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants have Boc and product has amine but no Boc
                reactant_has_boc = any(checker.check_fg("Boc", r) for r in reactants)
                product_has_boc = checker.check_fg("Boc", product)
                product_has_amine = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                print(
                    f"Manual check - Reactant has Boc: {reactant_has_boc}, Product has Boc: {product_has_boc}, Product has amine: {product_has_amine}"
                )

                # If reactant has Boc, product doesn't have Boc but has amine, it's likely a Boc deprotection
                if reactant_has_boc and not product_has_boc and product_has_amine:
                    is_boc_deprotection = True
                    print(f"Detected Boc deprotection by functional group analysis")

            if is_boc_deprotection:
                # Consider depths 0, 1, 2, and 3 as late-stage
                if depth <= 3:
                    boc_deprotection_at_late_stage = True
                    print(f"Late-stage Boc deprotection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {boc_deprotection_at_late_stage}")
    return boc_deprotection_at_late_stage
