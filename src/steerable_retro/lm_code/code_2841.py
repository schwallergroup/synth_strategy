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
    Detects an alcohol-to-bromide activation strategy for nucleophilic substitution.
    """
    # Track if we found alcohol to bromide conversion
    alcohol_to_bromide_found = False

    def dfs_traverse(node):
        nonlocal alcohol_to_bromide_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for specific alcohol to bromide reactions directly
                specific_reactions = [
                    "PBr3 and alcohol to alkyl bromide",
                    "Appel reaction",
                    "Wohl-Ziegler bromination benzyl primary",
                    "Wohl-Ziegler bromination benzyl secondary",
                    "Wohl-Ziegler bromination benzyl tertiary",
                    "Wohl-Ziegler bromination allyl primary",
                    "Wohl-Ziegler bromination allyl secondary",
                    "Wohl-Ziegler bromination allyl tertiary",
                    "Alcohol to bromide_Other",
                ]

                for rxn_type in specific_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found specific reaction {rxn_type}: {rsmi}")
                        alcohol_to_bromide_found = True
                        return

                # Check for forward direction (alcohol → bromide)
                alcohol_in_reactants = False
                alcohol_reactant = None
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                    ):
                        alcohol_in_reactants = True
                        alcohol_reactant = reactant
                        print(f"Found alcohol in reactant: {reactant}")
                        break

                bromide_in_product = False
                if (
                    checker.check_fg("Primary halide", product)
                    or checker.check_fg("Secondary halide", product)
                    or checker.check_fg("Tertiary halide", product)
                ):
                    # Verify it's specifically a bromide
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "Br":
                                bromide_in_product = True
                                print(f"Found bromide in product: {product}")
                                break

                # Check for retrosynthetic direction (bromide → alcohol)
                bromide_in_reactants = False
                bromide_reactant = None
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        for atom in reactant_mol.GetAtoms():
                            if atom.GetSymbol() == "Br":
                                bromide_in_reactants = True
                                bromide_reactant = reactant
                                print(f"Found bromide in reactant: {reactant}")
                                break
                        if bromide_in_reactants:
                            break

                alcohol_in_product = False
                if (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                ):
                    alcohol_in_product = True
                    print(f"Found alcohol in product: {product}")

                # Confirm transformation in either direction
                if (alcohol_in_reactants and bromide_in_product) or (
                    bromide_in_reactants and alcohol_in_product
                ):
                    print(f"Confirmed alcohol-bromide interconversion: {rsmi}")
                    alcohol_to_bromide_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return alcohol_to_bromide_found
