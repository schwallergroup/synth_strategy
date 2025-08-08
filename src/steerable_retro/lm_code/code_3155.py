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
    This function detects a synthetic strategy where a carboxylic acid is synthesized
    from a primary alcohol through functional group transformations.
    """
    # Track if we've found the strategy
    strategy_detected = False

    # Track molecules and reactions in the pathway
    alcohol_molecules = []
    aldehyde_molecules = []
    carboxylic_acid_molecules = []
    oxidation_reactions = []

    def dfs_traverse(node, depth=0, parent_smiles=None):
        nonlocal strategy_detected

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for primary alcohol
            if checker.check_fg("Primary alcohol", mol_smiles):
                print(f"Found primary alcohol at depth {depth}: {mol_smiles}")
                alcohol_molecules.append((mol_smiles, depth))

            # Check for aldehyde (potential intermediate)
            if checker.check_fg("Aldehyde", mol_smiles):
                print(f"Found aldehyde at depth {depth}: {mol_smiles}")
                aldehyde_molecules.append((mol_smiles, depth))

            # Check for carboxylic acid
            if checker.check_fg("Carboxylic acid", mol_smiles):
                print(f"Found carboxylic acid at depth {depth}: {mol_smiles}")
                carboxylic_acid_molecules.append((mol_smiles, depth))

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for oxidation reactions
            if (
                checker.check_reaction("Oxidation of alcohol to carboxylic acid", rxn_smiles)
                or checker.check_reaction("Oxidation of aldehyde to carboxylic acid", rxn_smiles)
                or checker.check_reaction("Oxidation of primary alcohol to aldehyde", rxn_smiles)
                or checker.check_reaction("Oxidation of alkene to carboxylic acid", rxn_smiles)
                or checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rxn_smiles)
            ):
                print(f"Found oxidation reaction at depth {depth}: {rxn_smiles}")
                oxidation_reactions.append((rxn_smiles, depth))

                # Extract reactants and product to verify transformation
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    # Check if reactant has alcohol and product has carboxylic acid
                    for reactant in reactants:
                        if checker.check_fg("Primary alcohol", reactant) and checker.check_fg(
                            "Carboxylic acid", product
                        ):
                            print(f"Confirmed alcohol to carboxylic acid transformation")
                            strategy_detected = True
                        # Also check for aldehyde intermediate
                        elif checker.check_fg("Aldehyde", reactant) and checker.check_fg(
                            "Carboxylic acid", product
                        ):
                            print(f"Confirmed aldehyde to carboxylic acid transformation")
                            # This is part of the strategy, but we need to find the alcohol to aldehyde part too
                        elif checker.check_fg("Primary alcohol", reactant) and checker.check_fg(
                            "Aldehyde", product
                        ):
                            print(f"Confirmed alcohol to aldehyde transformation")
                            # This is part of the strategy, but we need to find the aldehyde to acid part too
                except Exception as e:
                    print(f"Error extracting reactants and product: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, mol_smiles if node["type"] == "mol" else parent_smiles)

    # Start traversal from the root
    dfs_traverse(route)

    # If we didn't find direct evidence in reactions, check if we have both alcohols and carboxylic acids
    if not strategy_detected and alcohol_molecules and carboxylic_acid_molecules:
        # Check if there's a path from alcohol to carboxylic acid
        # Sort by depth to check progression (lower depth = later in synthesis, higher depth = earlier)
        alcohol_molecules.sort(key=lambda x: x[1])  # Lower depth first (later stage)
        carboxylic_acid_molecules.sort(
            key=lambda x: x[1], reverse=True
        )  # Higher depth first (earlier stage)

        # In retrosynthesis, the alcohol should be at a later stage (lower depth) than the acid
        if alcohol_molecules[0][1] < carboxylic_acid_molecules[0][1]:
            print(
                "Detected alcohol to carboxylic acid transformation strategy based on pathway analysis"
            )
            print(
                f"Alcohol at depth {alcohol_molecules[0][1]}, Carboxylic acid at depth {carboxylic_acid_molecules[0][1]}"
            )

            # Try to verify by comparing molecules directly
            alcohol_mol = Chem.MolFromSmiles(alcohol_molecules[0][0])
            acid_mol = Chem.MolFromSmiles(carboxylic_acid_molecules[0][0])

            if alcohol_mol and acid_mol:
                # Check if the molecules are similar enough to suggest a transformation
                alcohol_carbon_count = sum(
                    1 for atom in alcohol_mol.GetAtoms() if atom.GetSymbol() == "C"
                )
                acid_carbon_count = sum(
                    1 for atom in acid_mol.GetAtoms() if atom.GetSymbol() == "C"
                )

                # If carbon counts are similar (allowing for some difference due to protecting groups, etc.)
                carbon_diff = abs(alcohol_carbon_count - acid_carbon_count)
                if carbon_diff <= 3:  # Allow for some structural differences
                    print(
                        f"Molecule comparison supports transformation (carbon difference: {carbon_diff})"
                    )
                    strategy_detected = True

    # Check for two-step pathway through aldehyde
    if (
        not strategy_detected
        and alcohol_molecules
        and aldehyde_molecules
        and carboxylic_acid_molecules
    ):
        alcohol_molecules.sort(key=lambda x: x[1])  # Lower depth first (later stage)
        aldehyde_molecules.sort(key=lambda x: x[1])  # Sort aldehydes
        carboxylic_acid_molecules.sort(
            key=lambda x: x[1], reverse=True
        )  # Higher depth first (earlier stage)

        # Check if depths are in correct order: alcohol (later) -> aldehyde (middle) -> acid (earlier)
        if alcohol_molecules[0][1] < aldehyde_molecules[0][1] < carboxylic_acid_molecules[0][1]:
            print("Detected two-step alcohol → aldehyde → carboxylic acid pathway")
            print(
                f"Alcohol at depth {alcohol_molecules[0][1]}, Aldehyde at depth {aldehyde_molecules[0][1]}, Carboxylic acid at depth {carboxylic_acid_molecules[0][1]}"
            )
            strategy_detected = True

    # Final check: if we have both alcohols and carboxylic acids in the right order, that's enough evidence
    if not strategy_detected and alcohol_molecules and carboxylic_acid_molecules:
        alcohol_molecules.sort(key=lambda x: x[1])  # Lower depth first (later stage)
        carboxylic_acid_molecules.sort(
            key=lambda x: x[1], reverse=True
        )  # Higher depth first (earlier stage)

        if alcohol_molecules[0][1] < carboxylic_acid_molecules[0][1]:
            print(
                "Detected alcohol to carboxylic acid transformation based on molecule presence and depth ordering"
            )
            print(
                f"Alcohol at depth {alcohol_molecules[0][1]}, Carboxylic acid at depth {carboxylic_acid_molecules[0][1]}"
            )
            strategy_detected = True

    if not strategy_detected:
        print("Alcohol to carboxylic acid strategy not detected")

    return strategy_detected
