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
    Detects a synthetic strategy where a trifluoromethyl group is maintained
    throughout the synthesis.
    """
    # First, check if the target molecule has a trifluoromethyl group
    target_has_cf3 = False
    cf3_maintained = True

    # Find the target molecule (the root of the tree)
    if route["type"] == "mol" and not route.get("in_stock", False):
        target_mol_smiles = route["smiles"]
        target_has_cf3 = checker.check_fg("Trifluoro group", target_mol_smiles)
        print(f"Target molecule: {target_mol_smiles}")
        print(f"Target has CF3: {target_has_cf3}")

    # If target doesn't have CF3, no need to check maintenance
    if not target_has_cf3:
        return False

    # Function to identify reagents that don't need to have CF3
    def is_reagent(smiles):
        # Common reagents by exact match
        common_reagents = [
            "O=P(Cl)(Cl)Cl",  # POCl3
            "ClP(Cl)(Cl)(Cl)(Cl)Cl",  # PCl5
            "O=S(Cl)Cl",  # SOCl2
            "C(=O)(Cl)Cl",  # COCl2
            "BrP(Br)(Br)(Br)(Br)Br",  # PBr5
            "O=P(Br)(Br)Br",  # POBr3
            "Cl",
            "Br",
            "I",
            "F",
            "O",
            "N",
            "HCl",
            "HBr",
            "HI",
            "NH3",
        ]

        if any(reagent == smiles for reagent in common_reagents):
            return True

        # Check for small molecules that are likely reagents
        mol = Chem.MolFromSmiles(smiles)
        if mol and Chem.Descriptors.HeavyAtomCount(mol) <= 3:
            return True

        return False

    # Track CF3 maintenance through the synthesis
    def dfs_traverse(node, depth=0, is_main_path=True):
        nonlocal cf3_maintained

        # Skip if we already determined CF3 is not maintained
        if not cf3_maintained:
            return

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)
            is_stock = node.get("in_stock", False)

            print(
                f"Depth {depth}, Molecule: {mol_smiles}, Has CF3: {has_cf3}, In Stock: {is_stock}, Main Path: {is_main_path}"
            )

            # Only check molecules in the main synthetic path for CF3 maintenance
            if is_main_path and not is_stock and not has_cf3 and not is_reagent(mol_smiles):
                cf3_maintained = False
                print(f"CF3 not maintained in main path molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            # For reaction nodes, check if CF3 is introduced in this step
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)

                # If product has CF3, at least one reactant should have CF3
                if product_has_cf3:
                    # Filter out reagents
                    main_reactants = [r for r in reactants_smiles if not is_reagent(r)]

                    reactants_with_cf3 = any(
                        checker.check_fg("Trifluoro group", r) for r in main_reactants
                    )

                    if (
                        not reactants_with_cf3 and main_reactants
                    ):  # Only check if we have main reactants
                        print(f"CF3 introduced in reaction: {rsmi}")
                        cf3_maintained = False

                # Determine which children are in the main synthetic path
                # In retrosynthesis, the product is the current node we're examining
                for child in node.get("children", []):
                    if child["type"] == "mol":
                        # The product is in the main path, reactants are not
                        child_is_main_path = child["smiles"] == product_smiles
                        dfs_traverse(child, depth + 1, child_is_main_path)
                    else:
                        dfs_traverse(child, depth + 1, is_main_path)

                # Skip the regular traversal since we've handled it specially
                return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing for non-reaction nodes or if reaction handling failed
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, is_main_path)

    # Start traversal
    dfs_traverse(route)

    return cf3_maintained and target_has_cf3
