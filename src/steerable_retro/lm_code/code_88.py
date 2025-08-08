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
    This function detects if the synthesis route involves a functional group interconversion
    (specifically methoxy to bromo conversion).
    """
    # Track both direct and multi-step conversions
    methoxy_to_bromo_conversion = False
    methoxy_molecules = []  # Track molecules with methoxy groups
    bromo_molecules = []  # Track molecules with bromo groups

    def dfs_traverse(node, depth=0):
        nonlocal methoxy_to_bromo_conversion

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol is None:
                print(f"Warning: Could not parse molecule SMILES: {mol_smiles}")
                return

            # Check for methoxy groups
            if "OC" in mol_smiles and checker.check_fg("Ether", mol_smiles):
                # More specific check for methoxy pattern
                if mol.HasSubstructMatch(Chem.MolFromSmarts("CO[#6]")):
                    methoxy_molecules.append(mol_smiles)
                    print(f"Found molecule with methoxy group: {mol_smiles}")

            # Check for bromo groups
            if "Br" in mol_smiles and (
                checker.check_fg("Aromatic halide", mol_smiles)
                or checker.check_fg("Primary halide", mol_smiles)
                or checker.check_fg("Secondary halide", mol_smiles)
                or checker.check_fg("Tertiary halide", mol_smiles)
                or checker.check_fg("Alkenyl halide", mol_smiles)
            ):
                bromo_molecules.append(mol_smiles)
                print(f"Found molecule with bromo group: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this reaction converts methoxy to bromo
                methoxy_in_reactants = False
                bromo_in_product = False

                # Check reactants for methoxy groups
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and "OC" in reactant and checker.check_fg("Ether", reactant):
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("CO[#6]")):
                            methoxy_in_reactants = True
                            print(f"Found methoxy group in reactant: {reactant}")
                            break

                # Check product for bromo groups
                product_mol = Chem.MolFromSmiles(product_part)
                if (
                    product_mol
                    and "Br" in product_part
                    and (
                        checker.check_fg("Aromatic halide", product_part)
                        or checker.check_fg("Primary halide", product_part)
                        or checker.check_fg("Secondary halide", product_part)
                        or checker.check_fg("Tertiary halide", product_part)
                        or checker.check_fg("Alkenyl halide", product_part)
                    )
                ):
                    bromo_in_product = True
                    print(f"Found bromo group in product: {product_part}")

                # Check if this is a methoxy to bromo conversion reaction
                if methoxy_in_reactants and bromo_in_product:
                    # Check for relevant reactions
                    relevant_reactions = [
                        "Cleavage of methoxy ethers to alcohols",
                        "Alcohol to bromide",
                        "Wohl-Ziegler bromination benzyl primary",
                        "Alkyl bromides from alcohols",
                        "Cleavage of alkoxy ethers to alcohols",
                        "Alcohol to chloride_Other",
                        "Ether cleavage to primary alcohol",
                        "Primary alkyl halide to alcohol",
                        "Aromatic bromination",
                        "Bromination",
                    ]

                    for reaction_type in relevant_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Found methoxy to bromo conversion reaction: {reaction_type}")
                            print(f"Reaction SMILES: {rsmi}")
                            methoxy_to_bromo_conversion = True
                            return

                    # If no specific reaction type matched, check if the reaction involves
                    # both methoxy disappearance and bromo appearance
                    if (
                        "OCH3" in reactants_part
                        and "Br" in product_part
                        and "OCH3" not in product_part
                    ):
                        print(
                            "Detected potential methoxy to bromo conversion based on SMILES patterns"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        methoxy_to_bromo_conversion = True
                        return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start the traversal
    dfs_traverse(route)

    # Check for multi-step conversion by analyzing tracked molecules
    if not methoxy_to_bromo_conversion and methoxy_molecules and bromo_molecules:
        # If we have molecules with methoxy groups and molecules with bromo groups
        # in the synthesis route, this might indicate a multi-step conversion
        print(f"Detected potential multi-step methoxy to bromo conversion")
        print(f"Molecules with methoxy: {methoxy_molecules}")
        print(f"Molecules with bromo: {bromo_molecules}")
        methoxy_to_bromo_conversion = True

    print(
        f"Methoxy to bromo functional group interconversion detected: {methoxy_to_bromo_conversion}"
    )
    return methoxy_to_bromo_conversion
