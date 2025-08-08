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
    Detects if the final step (depth 0 or 1) involves Boc deprotection.
    """
    final_step_has_boc_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_has_boc_deprotection

        if node["type"] == "reaction" and (depth == 0 or depth == 1):
            # This is the final or near-final reaction step
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                reactants = reactants_str.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is any type of Boc deprotection reaction
                boc_deprotection_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                ]

                for rxn_type in boc_deprotection_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction")
                        final_step_has_boc_deprotection = True
                        return

                # Alternative check: look for Boc group in reactants but not in product
                boc_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Boc", reactant):
                        print(f"Found Boc group in reactant: {reactant}")
                        boc_in_reactants = True
                        break

                boc_in_product = checker.check_fg("Boc", product)
                if boc_in_product:
                    print(f"Found Boc group in product: {product}")

                if boc_in_reactants and not boc_in_product:
                    print("Detected Boc group removal (reactant has Boc, product doesn't)")

                    # Check if the reactant has a protected amine (N-Boc)
                    for reactant in reactants:
                        if checker.check_fg("Boc", reactant) and checker.check_fg(
                            "Primary amine", reactant
                        ):
                            print(
                                "Reactant has both Boc and amine groups - likely Boc deprotection"
                            )
                            final_step_has_boc_deprotection = True
                            return

                        # Check for secondary or tertiary amines too
                        if checker.check_fg("Boc", reactant) and (
                            checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                        ):
                            print(
                                "Reactant has both Boc and amine groups - likely Boc deprotection"
                            )
                            final_step_has_boc_deprotection = True
                            return

                # Additional check for any reaction that might be Boc deprotection
                if "Boc" in rsmi.lower() and any(
                    term in rsmi.lower() for term in ["deprotect", "cleav", "remov"]
                ):
                    print("Detected potential Boc deprotection based on SMILES string content")
                    final_step_has_boc_deprotection = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if (
                not final_step_has_boc_deprotection
            ):  # Stop traversal if we already found what we need
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {final_step_has_boc_deprotection}")
    return final_step_has_boc_deprotection
