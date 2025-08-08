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
    Detects a strategy involving cyclopropane formation followed by ring opening.
    """
    # Track molecules with cyclopropane and their synthetic history
    cyclopropane_molecules = set()
    ring_opening_reactions = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains a cyclopropane
            if checker.check_ring("cyclopropane", mol_smiles):
                print(f"Depth {depth}: Found molecule with cyclopropane: {mol_smiles}")
                cyclopropane_molecules.add(mol_smiles)

        elif node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # In retrosynthetic direction: product -> reactants
            # Check if product has cyclopropane but at least one reactant doesn't
            product_has_cyclopropane = checker.check_ring("cyclopropane", product_smiles)

            if product_has_cyclopropane:
                # Check if any reactant doesn't have cyclopropane (cyclopropane formation)
                for reactant in reactants_smiles:
                    if not checker.check_ring("cyclopropane", reactant):
                        print(f"Depth {depth}: Found cyclopropane formation reaction:")
                        print(f"  Reactants: {reactants_smiles}")
                        print(f"  Product: {product_smiles}")
                        # Add product to our tracked cyclopropane molecules
                        cyclopropane_molecules.add(product_smiles)
                        break

            # Check for ring opening: product has cyclopropane, reactant doesn't
            # In retrosynthetic direction, we're going from product to reactants
            if product_has_cyclopropane:
                for reactant in reactants_smiles:
                    if not checker.check_ring("cyclopropane", reactant):
                        print(f"Depth {depth}: Found cyclopropane ring opening reaction:")
                        print(f"  Product with ring: {product_smiles}")
                        print(f"  Reactant without ring: {reactant}")
                        ring_opening_reactions.add(rsmi)
                        break

            # Also check for functional group transformations on cyclopropane
            if product_has_cyclopropane:
                for reactant in reactants_smiles:
                    if checker.check_ring("cyclopropane", reactant):
                        # Check for ester hydrolysis (common cyclopropane transformation)
                        if (
                            checker.check_fg("Carboxylic acid", product_smiles)
                            and checker.check_fg("Ester", reactant)
                            and checker.check_reaction(
                                "Ester saponification (alkyl deprotection)", rsmi
                            )
                        ):
                            print(
                                f"Depth {depth}: Found cyclopropane functional group transformation:"
                            )
                            print(f"  Reactant (ester): {reactant}")
                            print(f"  Product (acid): {product_smiles}")
                            ring_opening_reactions.add(rsmi)
                            break

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both formation and opening
    has_formation = len(cyclopropane_molecules) > 0
    has_opening = len(ring_opening_reactions) > 0

    print(f"Found cyclopropane molecules: {len(cyclopropane_molecules)}")
    print(f"Found ring opening reactions: {len(ring_opening_reactions)}")

    return has_formation and has_opening
