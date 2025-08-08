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
    This function detects a synthetic strategy involving early-stage imine formation
    followed by late-stage aromatic C-O bond formations using organolithium reagents.
    """
    # Initialize tracking variables
    imine_formation_depth = -1  # -1 means not found
    organolithium_co_reactions = []  # Store depths of organolithium C-O bond formations
    max_depth = 0  # Track maximum depth to determine early vs late stage

    def dfs_traverse(node, depth=0):
        nonlocal imine_formation_depth, organolithium_co_reactions, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for imine formation
            if checker.check_fg("Substituted imine", product_smiles) and not any(
                checker.check_fg("Substituted imine", r) for r in reactants_smiles
            ):
                # Check if one reactant has aldehyde and another has amine
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)
                has_amine = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants_smiles
                )

                if has_aldehyde and has_amine:
                    print(f"Detected imine formation reaction at depth {depth}")
                    if imine_formation_depth == -1 or depth > imine_formation_depth:
                        imine_formation_depth = depth

            # Check for organolithium reagent and C-O bond formation
            has_organolithium = any(
                checker.check_fg("Alkyl lithium", r) or checker.check_fg("Aryl lithium", r)
                for r in reactants_smiles
            )

            if has_organolithium:
                print(f"Detected organolithium reagent at depth {depth}")

                # Check for specific organolithium reactions that form C-O bonds
                if checker.check_reaction("Directed ortho metalation of arenes", rsmi):
                    print(f"Detected C-O bond formation via ortho metalation at depth {depth}")
                    organolithium_co_reactions.append(depth)
                elif checker.check_reaction("Aryllithium cross-coupling", rsmi):
                    print(f"Detected aryllithium cross-coupling at depth {depth}")
                    organolithium_co_reactions.append(depth)
                elif checker.check_reaction("Carboxylic acid from Li and CO2", rsmi):
                    print(f"Detected organolithium carbonylation to acid at depth {depth}")
                    organolithium_co_reactions.append(depth)
                elif checker.check_reaction("Ketone from Li and CO2", rsmi):
                    print(f"Detected organolithium carbonylation to ketone at depth {depth}")
                    organolithium_co_reactions.append(depth)
                elif checker.check_reaction("Ketone from Li, Grignard and CO2", rsmi):
                    print(f"Detected organolithium carbonylation with Grignard at depth {depth}")
                    organolithium_co_reactions.append(depth)
                elif checker.check_reaction("Ketone from Li, halide and CO2", rsmi):
                    print(f"Detected organolithium carbonylation with halide at depth {depth}")
                    organolithium_co_reactions.append(depth)
                else:
                    # Check for general C-O bond formation
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    reactants_mols = [
                        Chem.MolFromSmiles(r)
                        for r in reactants_smiles
                        if not (
                            checker.check_fg("Alkyl lithium", r)
                            or checker.check_fg("Aryl lithium", r)
                        )
                    ]

                    # Check for phenol or ether formation
                    if (
                        checker.check_fg("Phenol", product_smiles)
                        or checker.check_fg("Ether", product_smiles)
                    ) and not all(
                        checker.check_fg("Phenol", r) or checker.check_fg("Ether", r)
                        for r in reactants_smiles
                    ):
                        print(f"Detected C-O bond formation at depth {depth}")
                        organolithium_co_reactions.append(depth)

                    # Check for methoxy group introduction
                    elif product_mol and reactants_mols and all(reactants_mols):
                        methoxy_pattern = Chem.MolFromSmarts("cOC")
                        if product_mol and methoxy_pattern:
                            product_methoxy_count = len(
                                product_mol.GetSubstructMatches(methoxy_pattern)
                            )

                            reactant_methoxy_count = 0
                            for r_mol in reactants_mols:
                                if r_mol:  # Ensure molecule is valid
                                    reactant_methoxy_count += len(
                                        r_mol.GetSubstructMatches(methoxy_pattern)
                                    )

                            if product_methoxy_count > reactant_methoxy_count:
                                print(f"Detected methoxy introduction at depth {depth}")
                                organolithium_co_reactions.append(depth)

            # Check for other reactions that might involve organolithium and C-O bond formation
            elif checker.check_reaction("O-methylation", rsmi) or checker.check_reaction(
                "Methylation of OH with DMS", rsmi
            ):
                # Check if there's evidence of organolithium involvement
                for reactant in reactants_smiles:
                    if (
                        checker.check_fg("Phenol", reactant)
                        or checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                    ):
                        print(
                            f"Detected potential organolithium-mediated O-methylation at depth {depth}"
                        )
                        organolithium_co_reactions.append(depth)
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    has_imine_formation = imine_formation_depth != -1
    has_organolithium_co = len(organolithium_co_reactions) >= 1

    # Check if imine formation is early stage (higher depth) and
    # organolithium reactions are late stage (lower depth)
    early_imine = has_imine_formation and imine_formation_depth >= max_depth / 2
    late_organolithium = has_organolithium_co and any(
        depth <= max_depth / 2 for depth in organolithium_co_reactions
    )

    strategy_present = early_imine and late_organolithium

    print(
        f"Strategy detection results: imine_formation={has_imine_formation} at depth {imine_formation_depth}, "
        f"organolithium_co_reactions={len(organolithium_co_reactions)} at depths {organolithium_co_reactions}, "
        f"max_depth={max_depth}, early_imine={early_imine}, late_organolithium={late_organolithium}"
    )

    return strategy_present
