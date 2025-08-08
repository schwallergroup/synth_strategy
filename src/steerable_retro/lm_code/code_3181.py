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
    This function detects a synthetic strategy involving disconnection of a tertiary amine
    to form a primary amine.
    """
    tertiary_amine_disconnection_found = False

    def dfs_traverse(node):
        nonlocal tertiary_amine_disconnection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a tertiary amine disconnection reaction
            if checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi):
                print("Found tertiary amine disconnection reaction via reaction type check")
                tertiary_amine_disconnection_found = True
            else:
                # If specific reaction check fails, check by functional group transformation
                reactants_str = rsmi.split(">")[0]
                products_str = rsmi.split(">")[-1]

                # In retrosynthesis, we're looking for tertiary amine in product being disconnected to primary amine in reactant
                reactants = reactants_str.split(".")
                products = products_str.split(".")

                # Check if any product has a tertiary amine and any reactant has a primary amine
                has_primary_amine_reactant = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )
                has_tertiary_amine_product = any(
                    checker.check_fg("Tertiary amine", p) for p in products
                )

                if has_primary_amine_reactant and has_tertiary_amine_product:
                    print("Found potential tertiary amine disconnection reaction")

                    # Try to verify this is a disconnection by checking atom mapping
                    try:
                        # Find the product with tertiary amine
                        tertiary_product = next(
                            p for p in products if checker.check_fg("Tertiary amine", p)
                        )
                        # Find the reactant with primary amine
                        primary_reactant = next(
                            r for r in reactants if checker.check_fg("Primary amine", r)
                        )

                        # Get the atom indices for the amines
                        tertiary_product_mol = Chem.MolFromSmiles(tertiary_product)
                        primary_reactant_mol = Chem.MolFromSmiles(primary_reactant)

                        if tertiary_product_mol and primary_reactant_mol:
                            # Check if this could be a hydrogenolysis or similar disconnection
                            # by verifying the atom count difference is reasonable
                            tertiary_atom_count = tertiary_product_mol.GetNumAtoms()
                            primary_atom_count = primary_reactant_mol.GetNumAtoms()

                            # Tertiary amine should have more atoms than primary amine
                            # and the difference should be reasonable for alkyl groups
                            if tertiary_atom_count > primary_atom_count:
                                print(
                                    f"Atom count check passed: tertiary ({tertiary_atom_count}) > primary ({primary_atom_count})"
                                )
                                tertiary_amine_disconnection_found = True
                    except (StopIteration, Exception) as e:
                        print(f"Error during atom mapping check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return tertiary_amine_disconnection_found
