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
    Detects a strategy where key functional groups (nitrile, halogens) are preserved
    while the molecular scaffold undergoes significant modifications.
    """
    # Track functional groups across the synthesis
    functional_groups = {
        "nitrile": {"present_at_start": False, "preserved": True},
        "halogen": {"present_at_start": False, "preserved": True},
    }

    # Track scaffold modifications
    scaffold_modified = False

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_modified

        if node["type"] == "mol":
            # Check for functional groups in the final product (depth 0)
            if depth == 0 and node["smiles"]:
                mol_smiles = node["smiles"]
                print(f"Checking final product: {mol_smiles}")

                # Check for nitrile in final product
                if checker.check_fg("Nitrile", mol_smiles):
                    functional_groups["nitrile"]["present_at_start"] = True
                    print("Nitrile found in final product")

                # Check for halogens in final product
                halogen_groups = [
                    "Aromatic halide",
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Alkenyl halide",
                    "Haloalkyne",
                ]
                for halogen_group in halogen_groups:
                    if checker.check_fg(halogen_group, mol_smiles):
                        functional_groups["halogen"]["present_at_start"] = True
                        print(f"{halogen_group} found in final product")
                        break

        elif node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for functional group preservation
                # Check nitrile preservation
                if any(r and checker.check_fg("Nitrile", r) for r in reactants_smiles):
                    if not checker.check_fg("Nitrile", product_smiles):
                        functional_groups["nitrile"]["preserved"] = False
                        print(f"Nitrile not preserved in reaction at depth {depth}")

                # Check halogen preservation
                halogen_groups = [
                    "Aromatic halide",
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Alkenyl halide",
                    "Haloalkyne",
                ]
                if any(
                    r and any(checker.check_fg(hg, r) for hg in halogen_groups)
                    for r in reactants_smiles
                ):
                    if not any(checker.check_fg(hg, product_smiles) for hg in halogen_groups):
                        functional_groups["halogen"]["preserved"] = False
                        print(f"Halogen not preserved in reaction at depth {depth}")

                # Check for scaffold modifications
                # Check for sulfone/sulfonamide modifications
                sulfur_groups = ["Sulfone", "Sulfonamide", "Sulfonate", "Sulfoxide"]
                for sg in sulfur_groups:
                    if any(r and checker.check_fg(sg, r) for r in reactants_smiles):
                        if not checker.check_fg(sg, product_smiles):
                            scaffold_modified = True
                            print(f"Detected scaffold modification at depth {depth}: {sg} removal")

                # Check for ring modifications
                for reactant_smiles in reactants_smiles:
                    if reactant_smiles and product_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if reactant_mol and product_mol:
                            reactant_ring_count = len(AllChem.GetSSSR(reactant_mol))
                            product_ring_count = len(AllChem.GetSSSR(product_mol))

                            if reactant_ring_count != product_ring_count:
                                scaffold_modified = True
                                print(
                                    f"Detected scaffold modification at depth {depth}: ring count change ({reactant_ring_count} to {product_ring_count})"
                                )
                                break

                # Check for reactions that typically modify scaffolds
                scaffold_modifying_reactions = [
                    "Diels-Alder",
                    "Retro-Diels-Alder from oxazole",
                    "Ring opening of epoxide with amine",
                    "Acetal hydrolysis to diol",
                    "Acetal hydrolysis to aldehyde",
                    "Ketal hydrolysis to ketone",
                    "Intramolecular transesterification/Lactone formation",
                    "Formation of NOS Heterocycles",
                ]

                for reaction in scaffold_modifying_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        scaffold_modified = True
                        print(f"Detected scaffold modification at depth {depth}: {reaction}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    nitrile_strategy = (
        functional_groups["nitrile"]["present_at_start"]
        and functional_groups["nitrile"]["preserved"]
    )

    halogen_strategy = (
        functional_groups["halogen"]["present_at_start"]
        and functional_groups["halogen"]["preserved"]
    )

    strategy_present = (nitrile_strategy or halogen_strategy) and scaffold_modified

    if strategy_present:
        print("Detected functional group preservation with scaffold modification strategy")
        if nitrile_strategy:
            print("- Nitrile group preserved throughout synthesis")
        if halogen_strategy:
            print("- Halogen group preserved throughout synthesis")
    else:
        if not (nitrile_strategy or halogen_strategy):
            print("No preserved functional groups detected")
        if not scaffold_modified:
            print("No significant scaffold modifications detected")

    return strategy_present
