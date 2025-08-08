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
    Detects if the synthesis involves a step where a halogen (F, Cl, Br, I)
    is removed from a molecule.
    """
    has_halogen_removal = False

    def dfs_traverse(node):
        nonlocal has_halogen_removal

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a known dehalogenation reaction
            if checker.check_reaction("Aromatic dehalogenation", rsmi) or checker.check_reaction(
                "Dehalogenation", rsmi
            ):
                print(f"Detected known dehalogenation reaction: {rsmi}")
                has_halogen_removal = True
                return

            # Manual check for halogen removal
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Count halogens in reactants
            reactants = reactants_part.split(".")
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            reactant_mols = [m for m in reactant_mols if m is not None]

            total_halogens_reactants = 0
            for mol in reactant_mols:
                for halogen in ["F", "Cl", "Br", "I"]:
                    pattern = Chem.MolFromSmarts(f"[{halogen}]")
                    if mol.HasSubstructMatch(pattern):
                        matches = mol.GetSubstructMatches(pattern)
                        total_halogens_reactants += len(matches)

            # Count halogens in products
            products = products_part.split(".")
            product_mols = [Chem.MolFromSmiles(p) for p in products if p]
            product_mols = [m for m in product_mols if m is not None]

            total_halogens_products = 0
            for mol in product_mols:
                for halogen in ["F", "Cl", "Br", "I"]:
                    pattern = Chem.MolFromSmarts(f"[{halogen}]")
                    if mol.HasSubstructMatch(pattern):
                        matches = mol.GetSubstructMatches(pattern)
                        total_halogens_products += len(matches)

            # Check if halogens were removed
            if total_halogens_reactants > total_halogens_products:
                print(
                    f"Detected halogen removal: {total_halogens_reactants} halogens in reactants, {total_halogens_products} in products"
                )
                print(f"Reaction: {rsmi}")
                has_halogen_removal = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_halogen_removal
