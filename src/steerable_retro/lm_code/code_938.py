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
    This function detects a linear synthesis strategy with multiple C-O bond
    manipulations while maintaining aromatic rings with specific substituents.
    """
    # Initialize tracking variables
    co_bond_manipulations = 0
    has_fluorinated_aromatic = False
    is_linear_synthesis = True

    def dfs_traverse(node, depth=0):
        nonlocal co_bond_manipulations, has_fluorinated_aromatic, is_linear_synthesis

        if node["type"] == "mol":
            # Check for fluorinated aromatic in molecules
            mol_smiles = node["smiles"]
            if checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles:
                # Verify F is connected to an aromatic system
                mol = Chem.MolFromSmiles(mol_smiles)
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F":
                        neighbor = atom.GetNeighbors()[0]  # F has only one neighbor
                        if neighbor.GetIsAromatic():
                            has_fluorinated_aromatic = True
                            print(f"Found fluorinated aromatic in: {mol_smiles}")
                            break

        elif node["type"] == "reaction":
            # Check for C-O bond manipulations
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction: {rsmi}")

                # Check if synthesis is linear (only one non-commercial reactant)
                non_commercial_children = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_commercial_children += 1

                if non_commercial_children > 1:
                    is_linear_synthesis = False
                    print(
                        f"Non-linear synthesis detected: {non_commercial_children} non-commercial reactants"
                    )

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for various C-O bond manipulations by reaction type
                if (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Alcohol to ether", rsmi)
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                    or checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                ):
                    co_bond_manipulations += 1
                    print(f"C-O bond manipulation detected by reaction type in: {rsmi}")

                # Check for C-O bond manipulations by examining functional group changes
                else:
                    reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    reactant_has_ether = any(checker.check_fg("Ether", r) for r in reactants)
                    reactant_has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )
                    reactant_has_carboxylic = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )

                    product_has_ester = checker.check_fg("Ester", product)
                    product_has_ether = checker.check_fg("Ether", product)
                    product_has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )
                    product_has_carboxylic = checker.check_fg("Carboxylic acid", product)

                    # Check for lactone formation (intramolecular esterification)
                    reactant_has_lactone = any(
                        checker.check_ring("oxolane", r) or checker.check_ring("oxane", r)
                        for r in reactants
                    )
                    product_has_lactone = checker.check_ring(
                        "oxolane", product
                    ) or checker.check_ring("oxane", product)

                    if (
                        (reactant_has_ester and not product_has_ester)
                        or (not reactant_has_ester and product_has_ester)
                        or (reactant_has_ether and not product_has_ether)
                        or (not reactant_has_ether and product_has_ether)
                        or (reactant_has_alcohol and not product_has_alcohol)
                        or (not reactant_has_alcohol and product_has_alcohol)
                        or (reactant_has_carboxylic and not product_has_carboxylic)
                        or (not reactant_has_carboxylic and product_has_carboxylic)
                        or (reactant_has_lactone and not product_has_lactone)
                        or (not reactant_has_lactone and product_has_lactone)
                    ):
                        co_bond_manipulations += 1
                        print(f"C-O bond manipulation detected through FG change in: {rsmi}")

                    # Check for intramolecular cyclization reactions that form C-O bonds
                    elif checker.check_reaction(
                        "Intramolecular transesterification/Lactone formation", rsmi
                    ):
                        co_bond_manipulations += 1
                        print(
                            f"C-O bond manipulation detected through lactone formation in: {rsmi}"
                        )

                    # Check for reactions that might involve C-O bond formation/cleavage
                    elif (
                        checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                        or checker.check_reaction(
                            "Reduction of aldehydes and ketones to alcohols", rsmi
                        )
                        or checker.check_reaction(
                            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                            rsmi,
                        )
                    ):
                        co_bond_manipulations += 1
                        print(f"C-O bond manipulation detected through related reaction in: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if strategy is present - relaxed to 2+ C-O manipulations
    strategy_present = (
        co_bond_manipulations >= 2 and has_fluorinated_aromatic and is_linear_synthesis
    )

    print(f"C-O bond manipulations: {co_bond_manipulations}")
    print(f"Has fluorinated aromatic: {has_fluorinated_aromatic}")
    print(f"Is linear synthesis: {is_linear_synthesis}")
    print(f"Strategy present: {strategy_present}")

    return strategy_present
