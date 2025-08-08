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
    This function detects a synthetic strategy involving incorporation of a
    dichlorophenyl group via an isocyanate intermediate.
    """
    # Track if we found the strategy
    strategy_found = False

    # Patterns for dichlorophenyl groups (1,2-, 1,3-, and 1,4-dichlorophenyl)
    dichlorophenyl_patterns = ["c1c(Cl)cc(Cl)cc1", "c1c(Cl)ccc(Cl)c1", "c1cc(Cl)c(Cl)cc1"]

    def dfs_traverse(node, depth=0):
        nonlocal strategy_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            try:
                # Check for dichlorophenyl-containing isocyanate in reactants
                dichlorophenyl_isocyanate_in_reactants = False
                for reactant_smiles in reactants_smiles:
                    # Check if this reactant has dichlorophenyl pattern
                    has_dichlorophenyl = any(
                        pattern in reactant_smiles for pattern in dichlorophenyl_patterns
                    ) or (
                        checker.check_fg("Aromatic halide", reactant_smiles)
                        and reactant_smiles.count("Cl") >= 2
                    )

                    # Check if this reactant has isocyanate
                    has_isocyanate = checker.check_fg("Isocyanate", reactant_smiles)

                    print(f"Reactant: {reactant_smiles}")
                    print(f"  Has dichlorophenyl: {has_dichlorophenyl}")
                    print(f"  Has isocyanate: {has_isocyanate}")

                    if has_dichlorophenyl and has_isocyanate:
                        print(
                            f"Found dichlorophenyl-containing isocyanate in reactant: {reactant_smiles}"
                        )
                        dichlorophenyl_isocyanate_in_reactants = True
                        break

                # Check if the product contains a structure derived from dichlorophenyl-isocyanate
                if dichlorophenyl_isocyanate_in_reactants:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        has_dichlorophenyl_in_product = any(
                            pattern in product_smiles for pattern in dichlorophenyl_patterns
                        ) or (
                            checker.check_fg("Aromatic halide", product_smiles)
                            and product_smiles.count("Cl") >= 2
                        )

                        print(f"Product: {product_smiles}")
                        print(f"  Has dichlorophenyl: {has_dichlorophenyl_in_product}")

                        # Check for common products of isocyanate reactions
                        has_urea = checker.check_fg("Urea", product_smiles)
                        has_carbamate = checker.check_fg("Carbamic ester", product_smiles)
                        has_amide = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )

                        print(f"  Has urea: {has_urea}")
                        print(f"  Has carbamate: {has_carbamate}")
                        print(f"  Has amide: {has_amide}")

                        if has_dichlorophenyl_in_product and (
                            has_urea or has_carbamate or has_amide
                        ):
                            print(
                                f"Found incorporation of dichlorophenyl group via isocyanate in product: {product_smiles}"
                            )
                            strategy_found = True

                        # Check if the reaction is a known isocyanate reaction
                        isocyanate_reaction = (
                            checker.check_reaction(
                                "Urea synthesis via isocyanate and primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Urea synthesis via isocyanate and secondary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Urea synthesis via isocyanate and diazo", rsmi
                            )
                            or checker.check_reaction(
                                "Urea synthesis via isocyanate and sulfonamide", rsmi
                            )
                            or checker.check_reaction("Carbamic ester", rsmi)
                        )

                        print(f"  Is known isocyanate reaction: {isocyanate_reaction}")

                        if has_dichlorophenyl_in_product and isocyanate_reaction:
                            print(
                                f"Found dichlorophenyl incorporation via isocyanate reaction: {rsmi}"
                            )
                            strategy_found = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Strategy found: {strategy_found}")
    return strategy_found
