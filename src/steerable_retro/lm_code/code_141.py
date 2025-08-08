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
    This function detects the use of boronic acid as a key intermediate
    for subsequent coupling reactions.
    """
    boronic_acid_formed = False
    boronic_acid_used = False

    def dfs_traverse(node, depth=0):
        nonlocal boronic_acid_formed, boronic_acid_used

        # Process reaction nodes
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid formation reactions
                # In retrosynthesis, this means the product contains boronic acid
                product_mol = product
                if checker.check_fg("Boronic acid", product_mol):
                    print(f"Depth {depth}: Found boronic acid in product: {product_mol}")

                    # Check if this is a formation reaction (not just a modification of existing boronic acid)
                    has_boronic_acid_in_reactants = any(
                        checker.check_fg("Boronic acid", r) for r in reactants
                    )
                    has_boronic_ester_in_reactants = any(
                        checker.check_fg("Boronic ester", r) for r in reactants
                    )

                    if not has_boronic_acid_in_reactants and not has_boronic_ester_in_reactants:
                        boronic_acid_formed = True
                        print(f"Depth {depth}: Confirmed boronic acid formation")
                    elif has_boronic_ester_in_reactants:
                        # Conversion from boronic ester to boronic acid also counts as formation
                        boronic_acid_formed = True
                        print(f"Depth {depth}: Confirmed boronic acid formation from ester")

                # Check for boronic acid usage in coupling reactions
                # In retrosynthesis, this means a reactant contains boronic acid
                for reactant in reactants:
                    if checker.check_fg("Boronic acid", reactant):
                        print(f"Depth {depth}: Found boronic acid in reactant: {reactant}")

                        # Check if this is a coupling reaction
                        if (
                            checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                            or checker.check_reaction(
                                "Suzuki coupling with boronic acids OTf", rsmi
                            )
                            or checker.check_reaction("Chan-Lam alcohol", rsmi)
                            or checker.check_reaction("Chan-Lam amine", rsmi)
                            or checker.check_reaction(
                                "Petasis reaction with amines and boronic acids", rsmi
                            )
                        ):
                            boronic_acid_used = True
                            print(
                                f"Depth {depth}: Confirmed boronic acid used in coupling reaction"
                            )
                            break

                        # If specific reaction check fails, look for general coupling patterns
                        if not boronic_acid_used:
                            # Check if product doesn't have boronic acid (consumed in reaction)
                            if not checker.check_fg("Boronic acid", product):
                                boronic_acid_used = True
                                print(
                                    f"Depth {depth}: Boronic acid consumed in reaction (likely coupling)"
                                )

        # Process children nodes
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(
        f"Final result - Boronic acid formed: {boronic_acid_formed}, Boronic acid used: {boronic_acid_used}"
    )
    return boronic_acid_formed and boronic_acid_used
