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
    Detects if the synthesis introduces an alkyl chain with terminal functionality
    (like an ester group) to an aromatic core.
    """
    chain_introduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal chain_introduction_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an ether synthesis or similar reaction
                is_ether_synthesis = (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Chan-Lam etherification", rsmi)
                )

                # Check for aromatic core with ether linkage to a chain with terminal functionality
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    has_aromatic = any(atom.GetIsAromatic() for atom in product_mol.GetAtoms())
                    has_ether = checker.check_fg("Ether", product)
                    has_ester = checker.check_fg("Ester", product)
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
                    has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    )

                    print(
                        f"Product analysis - Aromatic: {has_aromatic}, Ether: {has_ether}, Ester: {has_ester}, Acid: {has_carboxylic_acid}, Alcohol: {has_alcohol}"
                    )

                    # Check if product has both aromatic core and terminal functionality connected by ether
                    if (
                        has_aromatic
                        and has_ether
                        and (has_ester or has_carboxylic_acid or has_alcohol)
                    ):
                        # Check if reactants don't have this complete pattern (it's being introduced)
                        reactants_with_pattern = 0
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if (
                                reactant_mol
                                and any(atom.GetIsAromatic() for atom in reactant_mol.GetAtoms())
                                and checker.check_fg("Ether", reactant)
                                and (
                                    checker.check_fg("Ester", reactant)
                                    or checker.check_fg("Carboxylic acid", reactant)
                                    or checker.check_fg("Primary alcohol", reactant)
                                    or checker.check_fg("Secondary alcohol", reactant)
                                    or checker.check_fg("Tertiary alcohol", reactant)
                                )
                            ):
                                reactants_with_pattern += 1
                                print(f"Reactant already has the pattern: {reactant}")

                        if reactants_with_pattern == 0:
                            # If we're introducing the pattern through an ether synthesis or at a reasonable stage
                            if (
                                is_ether_synthesis or depth <= 5
                            ):  # Allow for mid-stage modifications too
                                chain_introduction_found = True
                                print(
                                    f"Alkyl chain with terminal functionality introduction detected at depth {depth}"
                                )

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return chain_introduction_found
