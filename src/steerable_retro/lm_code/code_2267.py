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
    This function detects if the synthesis involves a late-stage deprotection
    of an allyl-protected alcohol.
    """
    # Track if we found the pattern
    found_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a shallow reaction (late in synthesis, depth 0 or 1)
            if depth <= 1:
                print(f"Checking late-stage reaction at depth {depth}: {rsmi}")

                # Check if any reactant contains allyl ether
                reactant_with_allyl_ether = None
                for reactant in reactants:
                    if checker.check_fg("Allyl", reactant) and checker.check_fg("Ether", reactant):
                        reactant_with_allyl_ether = reactant
                        print(f"Found reactant with allyl ether: {reactant}")
                        break

                # Check if product contains primary alcohol
                product_has_alcohol = checker.check_fg("Primary alcohol", product)
                if product_has_alcohol:
                    print(f"Product has primary alcohol: {product}")

                # Check if this is a deprotection reaction
                if reactant_with_allyl_ether and product_has_alcohol:
                    # Check if this is a deprotection reaction type
                    if checker.check_reaction(
                        "Hydroxyl benzyl deprotection", rsmi
                    ) or checker.check_reaction("Ether cleavage to primary alcohol", rsmi):
                        print("Found allyl ether deprotection in late-stage reaction")
                        found_deprotection = True
                    else:
                        # Additional check to confirm it's a deprotection
                        # Compare atom counts to ensure we're removing the allyl group
                        reactant_mol = Chem.MolFromSmiles(reactant_with_allyl_ether)
                        product_mol = Chem.MolFromSmiles(product)

                        if reactant_mol and product_mol:
                            # Allyl group typically has 3 carbons
                            if reactant_mol.GetNumAtoms() > product_mol.GetNumAtoms() + 2:
                                print("Atom count suggests deprotection (removing allyl group)")
                                found_deprotection = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_deprotection}")
    return found_deprotection
