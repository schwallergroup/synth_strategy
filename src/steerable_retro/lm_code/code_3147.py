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
    Detects if the synthesis uses an azide reduction followed by a cyclization step.
    """
    # Track if we found the key reactions
    found_azide_reduction = False
    found_cyclization = False
    azide_depth = -1
    cyclization_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_azide_reduction, found_cyclization, azide_depth, cyclization_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for azide reduction (N3 → NH2)
                if checker.check_reaction(
                    "Azide to amine reduction (Staudinger)", rsmi
                ) or checker.check_reaction("Azide to amine reduction", rsmi):
                    found_azide_reduction = True
                    azide_depth = depth
                    print(f"Found azide reduction at depth {depth}, rsmi: {rsmi}")
                elif any(
                    checker.check_fg("Azide", r) for r in reactants_smiles
                ) and checker.check_fg("Primary amine", product_smiles):
                    # Fallback check for azide reduction if specific reaction check fails
                    # Verify that azide is actually being reduced (not just present)
                    reactant_azide_count = sum(
                        1 for r in reactants_smiles if checker.check_fg("Azide", r)
                    )
                    product_azide_count = 1 if checker.check_fg("Azide", product_smiles) else 0

                    if reactant_azide_count > product_azide_count:
                        found_azide_reduction = True
                        azide_depth = depth
                        print(f"Found azide reduction (FG check) at depth {depth}, rsmi: {rsmi}")

                # Check for cyclization (formation of a new ring)
                # Convert to RDKit molecules for ring counting
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if all(r for r in reactants) and product:
                    # Count rings in reactants and product
                    reactant_rings = sum(len(Chem.GetSSSR(r)) for r in reactants if r)
                    product_rings = len(Chem.GetSSSR(product))

                    # Check for ring formation
                    if product_rings > reactant_rings:
                        # Check for specific cyclization reactions
                        cyclization_reactions = [
                            "Formation of NOS Heterocycles",
                            "Intramolecular amination (heterocycle formation)",
                            "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                            "Paal-Knorr pyrrole synthesis",
                            "Benzimidazole formation",
                            "Benzothiazole formation",
                            "Benzoxazole formation",
                        ]

                        is_specific_cyclization = any(
                            checker.check_reaction(rxn, rsmi) for rxn in cyclization_reactions
                        )

                        # Accept either specific reaction type or general ring formation
                        if is_specific_cyclization or product_rings > reactant_rings:
                            found_cyclization = True
                            cyclization_depth = depth
                            print(f"Found cyclization at depth {depth}, rsmi: {rsmi}")
                            print(f"  Rings: reactants={reactant_rings}, product={product_rings}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present: azide reduction followed by cyclization
    # In retrosynthetic traversal, higher depth means earlier stage in synthesis
    # So for azide reduction to happen before cyclization in forward synthesis,
    # azide_depth must be greater than cyclization_depth in retrosynthesis
    strategy_present = (
        found_azide_reduction and found_cyclization and azide_depth > cyclization_depth
    )

    if strategy_present:
        print(
            f"Detected azide reduction before cyclization strategy: azide_depth={azide_depth}, cyclization_depth={cyclization_depth}"
        )
    else:
        if found_azide_reduction and found_cyclization:
            print(
                f"Found both reactions but in wrong order: azide_depth={azide_depth}, cyclization_depth={cyclization_depth}"
            )
        elif not found_azide_reduction:
            print("No azide reduction found")
        elif not found_cyclization:
            print("No cyclization found")

    return strategy_present
