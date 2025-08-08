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
    This function detects a strategy involving early fluorination of an intermediate.
    """
    early_fluorination = False

    def dfs_traverse(node, depth=0):
        nonlocal early_fluorination

        if node["type"] == "reaction" and depth >= 1:  # Early in synthesis (depth >= 1)
            # Get reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if fluorine is introduced
                reactants_have_f = False
                for r_smi in reactants_smiles:
                    mol = Chem.MolFromSmiles(r_smi)
                    if mol and any(atom.GetSymbol() == "F" for atom in mol.GetAtoms()):
                        reactants_have_f = True
                        print(f"Reactant contains fluorine: {r_smi}")
                        break

                product_mol = Chem.MolFromSmiles(product_smiles)
                product_has_f = False
                if product_mol:
                    product_has_f = any(atom.GetSymbol() == "F" for atom in product_mol.GetAtoms())
                    if product_has_f:
                        print(f"Product contains fluorine: {product_smiles}")

                # Check if this is a fluorination reaction
                is_fluorination = checker.check_reaction("Fluorination", rsmi)
                is_aromatic_fluorination = checker.check_reaction("Aromatic fluorination", rsmi)

                # Check for C-F bond formation
                c_f_bond_formed = product_has_f and not reactants_have_f

                # Check if fluorine is transferred from reagent to product
                f_transfer = False
                if product_has_f and reactants_have_f:
                    # Check if the fluorine is part of a reagent and gets incorporated into the product
                    for r_smi in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol and any(atom.GetSymbol() == "F" for atom in r_mol.GetAtoms()):
                            # If this reactant is small and contains F, it might be a fluorinating agent
                            if r_mol.GetNumAtoms() < 15:  # Arbitrary threshold for "small molecule"
                                f_transfer = True
                                print(f"Potential fluorine transfer from reagent: {r_smi}")

                # Detect early fluorination
                if (c_f_bond_formed or f_transfer) and (
                    is_fluorination or is_aromatic_fluorination
                ):
                    print(f"Detected early fluorination at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    early_fluorination = True

                # Additional check for fluorination functional groups
                if not early_fluorination and product_has_f:
                    # Check if product contains fluorinated functional groups not present in reactants
                    has_trifluoro = checker.check_fg("Trifluoro group", product_smiles)
                    has_aromatic_halide = checker.check_fg("Aromatic halide", product_smiles)

                    reactants_have_trifluoro = any(
                        checker.check_fg("Trifluoro group", r) for r in reactants_smiles
                    )
                    reactants_have_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    )

                    if (has_trifluoro and not reactants_have_trifluoro) or (
                        has_aromatic_halide and not reactants_have_aromatic_halide
                    ):
                        print(f"Detected early fluorination via functional group at depth {depth}")
                        early_fluorination = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to detect early fluorination strategy")
    dfs_traverse(route)
    print(f"Early fluorination detected: {early_fluorination}")
    return early_fluorination
