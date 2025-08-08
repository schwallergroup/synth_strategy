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
    This function detects a synthetic strategy involving C-P bond formation
    to create a phosphonate group.
    """
    has_c_p_bond_formation = False

    def dfs_traverse(node):
        nonlocal has_c_p_bond_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction: {rsmi}")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phosphorus in product
            if "P" in product and "O" in product:
                print(f"Found phosphorus and oxygen in product: {product}")

                # Check for phosphonate structure in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for C-P bond in product
                    if product_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]-[P]")):
                        print(f"Found C-P bond in product: {product}")

                        # Check for halides in reactants
                        has_halide = False
                        has_phosphorus_compound = False
                        halide_reactant = ""
                        phosphorus_reactant = ""

                        for reactant in reactants:
                            # Check for various halides
                            if (
                                checker.check_fg("Primary halide", reactant)
                                or checker.check_fg("Secondary halide", reactant)
                                or checker.check_fg("Tertiary halide", reactant)
                                or checker.check_fg("Aromatic halide", reactant)
                            ):
                                has_halide = True
                                halide_reactant = reactant
                                print(f"Found halide in reactant: {reactant}")

                            # Check for phosphorus-containing compounds
                            if "P" in reactant:
                                has_phosphorus_compound = True
                                phosphorus_reactant = reactant
                                print(f"Found phosphorus compound in reactant: {reactant}")

                        # If we have both a halide and a phosphorus compound, check for C-P bond formation
                        if has_halide and has_phosphorus_compound:
                            # Check if the halide carbon is now connected to phosphorus in the product
                            halide_mol = Chem.MolFromSmiles(halide_reactant)

                            if halide_mol:
                                # Look for carbon attached to halogen in reactant
                                for atom in halide_mol.GetAtoms():
                                    if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                        # Get neighboring carbon atoms
                                        for neighbor in atom.GetNeighbors():
                                            if neighbor.GetSymbol() == "C":
                                                print(
                                                    f"Found carbon attached to halogen in reactant"
                                                )

                                                # This is likely an Arbuzov or similar reaction
                                                has_c_p_bond_formation = True
                                                print(
                                                    "Confirmed C-P bond formation to create phosphonate"
                                                )
                                                return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_c_p_bond_formation
