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
    This function detects a synthetic strategy involving late-stage introduction
    of a trifluoromethyl-containing sulfonamide group.
    """
    has_late_stage_cf3_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cf3_sulfonamide

        # Check reaction nodes at late stage (final or penultimate step)
        if node["type"] == "reaction" and depth <= 1:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    # Properly split reaction SMILES
                    parts = rsmi.split(">")
                    reactants_str = parts[0]
                    product_str = parts[-1]

                    # Skip if product is empty
                    if not product_str:
                        print("Empty product string, skipping")
                        return

                    reactants = reactants_str.split(".")
                    product_mol = Chem.MolFromSmiles(product_str)

                    # Check if this is a sulfonamide formation reaction
                    is_sulfonamide_reaction = (
                        checker.check_reaction(
                            "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                        )
                        or
                        # Check for sulfonyl chloride in reactants and S(=O)(=O)N pattern in product
                        (
                            any(checker.check_fg("Sulfonyl halide", r) for r in reactants if r)
                            and (
                                product_mol is not None
                                and product_mol.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)N"))
                            )
                        )
                    )
                    print(f"Is sulfonamide reaction: {is_sulfonamide_reaction}")

                    # Check for trifluoromethyl group in reactants
                    has_cf3_reactant = any(
                        checker.check_fg("Trifluoro group", r) for r in reactants if r
                    )
                    print(f"Has CF3 in reactants: {has_cf3_reactant}")

                    # Check if product contains both trifluoromethyl group and sulfonamide
                    has_cf3_product = checker.check_fg("Trifluoro group", product_str)

                    # Enhanced sulfonamide detection
                    has_sulfonamide_product = checker.check_fg("Sulfonamide", product_str)
                    if not has_sulfonamide_product and product_mol is not None:
                        has_sulfonamide_product = product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("S(=O)(=O)N")
                        )
                    print(f"Has CF3 in product: {has_cf3_product}")
                    print(f"Has sulfonamide in product: {has_sulfonamide_product}")

                    # Direct check for trifluoromethyl sulfonamide pattern
                    has_cf3_sulfonamide = False
                    if product_mol is not None:
                        # Check for CF3 group near sulfonamide
                        has_cf3_sulfonamide = product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("S(=O)(=O)N.*C(F)(F)F")
                        )
                    print(f"Direct check for CF3 sulfonamide: {has_cf3_sulfonamide}")

                    # Determine if this is a late-stage trifluoromethyl sulfonamide formation
                    if (
                        is_sulfonamide_reaction
                        and has_cf3_reactant
                        and (has_cf3_product and has_sulfonamide_product or has_cf3_sulfonamide)
                    ):
                        print(
                            f"Found late-stage trifluoromethyl sulfonamide formation at depth {depth}"
                        )
                        has_late_stage_cf3_sulfonamide = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_cf3_sulfonamide
