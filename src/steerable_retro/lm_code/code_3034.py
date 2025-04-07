#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a strategy where heterocyclic ring formation is followed by
    sequential functional group modifications
    """
    # Track ring formation and subsequent functional group modifications
    ring_formation_found = False
    ring_formation_depth = -1
    fg_modifications = []

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_found, ring_formation_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol is None or any(mol is None for mol in reactant_mols):
                    print(f"Warning: Could not parse some molecules in reaction: {rsmi}")
                    return

                # Count rings in reactants and product
                reactant_ring_count = sum(
                    mol.GetRingInfo().NumRings() for mol in reactant_mols if mol is not None
                )
                product_ring_count = (
                    product_mol.GetRingInfo().NumRings() if product_mol is not None else 0
                )

                print(
                    f"Depth: {depth}, Reactant rings: {reactant_ring_count}, Product rings: {product_ring_count}"
                )

                # Check for ring formation
                if product_ring_count > reactant_ring_count:
                    # Check if it's a heterocyclic ring formation
                    heterocyclic_formed = False

                    # Check for common heterocyclic rings
                    heterocyclic_rings = [
                        "furan",
                        "pyran",
                        "dioxane",
                        "tetrahydrofuran",
                        "pyrrole",
                        "pyridine",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "triazole",
                        "tetrazole",
                        "morpholine",
                        "thiomorpholine",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "thiophene",
                        "isoxazole",
                        "isothiazole",
                        "oxadiazole",
                        "thiadiazole",
                        "benzoxazole",
                        "benzothiazole",
                        "benzimidazole",
                    ]

                    for ring_name in heterocyclic_rings:
                        if checker.check_ring(ring_name, product_smiles) and not any(
                            checker.check_ring(ring_name, r) for r in reactants_smiles
                        ):
                            print(
                                f"Heterocyclic ring formation detected: {ring_name} at depth {depth}"
                            )
                            heterocyclic_formed = True
                            break

                    # If no specific ring matched, check for general heterocyclic formation
                    if not heterocyclic_formed:
                        for ring in product_mol.GetSSSR():
                            ring_atoms = set(
                                product_mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring
                            )
                            # Check for N, O, S atoms (atomic numbers 7, 8, 16)
                            if 7 in ring_atoms or 8 in ring_atoms or 16 in ring_atoms:
                                print(
                                    f"Generic heterocyclic ring formation detected at depth {depth}"
                                )
                                heterocyclic_formed = True
                                break

                    if heterocyclic_formed:
                        ring_formation_found = True
                        ring_formation_depth = depth

                # Check for functional group modifications using checker functions
                fg_mod_detected = False

                # Check for esterification
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    fg_modifications.append(("esterification", depth))
                    fg_mod_detected = True
                    print(f"Esterification detected at depth {depth}")

                # Check for reduction to alcohol
                if (
                    checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                ):
                    fg_modifications.append(("reduction", depth))
                    fg_mod_detected = True
                    print(f"Reduction detected at depth {depth}")

                # Check for mesylation
                if checker.check_fg("Mesylate", product_smiles) and not any(
                    checker.check_fg("Mesylate", r) for r in reactants_smiles
                ):
                    fg_modifications.append(("mesylation", depth))
                    fg_mod_detected = True
                    print(f"Mesylation detected at depth {depth}")

                # Check for other common functional group modifications
                if checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                ):
                    fg_modifications.append(("acylation", depth))
                    fg_mod_detected = True
                    print(f"Acylation detected at depth {depth}")

                if checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                ):
                    fg_modifications.append(("oxidation", depth))
                    fg_mod_detected = True
                    print(f"Oxidation detected at depth {depth}")

                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                    fg_modifications.append(("etherification", depth))
                    fg_mod_detected = True
                    print(f"Etherification detected at depth {depth}")

                if checker.check_reaction("Alkylation of amines", rsmi):
                    fg_modifications.append(("alkylation", depth))
                    fg_mod_detected = True
                    print(f"Alkylation detected at depth {depth}")

                # Check for protection/deprotection reactions
                if "protection" in rsmi.lower() or "deprotection" in rsmi.lower():
                    fg_modifications.append(("protection/deprotection", depth))
                    fg_mod_detected = True
                    print(f"Protection/deprotection detected at depth {depth}")

                # If no specific reaction matched but functional groups changed
                if not fg_mod_detected:
                    # Check for appearance/disappearance of common functional groups
                    common_fgs = [
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Carboxylic acid",
                        "Ester",
                        "Amide",
                        "Nitrile",
                        "Aldehyde",
                        "Ketone",
                    ]

                    for fg in common_fgs:
                        if (
                            checker.check_fg(fg, product_smiles)
                            and not any(checker.check_fg(fg, r) for r in reactants_smiles)
                        ) or (
                            any(checker.check_fg(fg, r) for r in reactants_smiles)
                            and not checker.check_fg(fg, product_smiles)
                        ):
                            fg_modifications.append((f"{fg} modification", depth))
                            print(f"{fg} modification detected at depth {depth}")
                            break

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Filter out FG modifications that occur at the same depth as ring formation
    # or at a higher depth (earlier stage in retrosynthesis)
    valid_fg_mods = [(mod, d) for mod, d in fg_modifications if d < ring_formation_depth]

    # Sort by depth to check if they're sequential (lower depth = later stage in retrosynthesis)
    valid_fg_mods.sort(key=lambda x: x[1])

    print(f"Valid FG modifications (after ring formation): {valid_fg_mods}")

    # Check if we have ring formation followed by at least 2 functional group modifications
    strategy_detected = ring_formation_found and len(valid_fg_mods) >= 2

    # Check if the modifications are sequential (consecutive depths or close)
    if strategy_detected and len(valid_fg_mods) >= 2:
        # Check if the modifications are reasonably sequential
        depths = [d for _, d in valid_fg_mods]
        is_sequential = True

        # Check if there are at least 2 modifications with depths less than ring_formation_depth
        # and these modifications are reasonably close to each other
        for i in range(len(depths) - 1):
            if depths[i + 1] - depths[i] > 3:  # Allow some flexibility in "sequential" definition
                is_sequential = False
                break

        strategy_detected = is_sequential

    print(f"Heterocyclic ring formation followed by FG modifications: {strategy_detected}")
    print(f"Ring formation depth: {ring_formation_depth if ring_formation_found else 'Not found'}")
    print(f"FG modifications: {fg_modifications}")
    print(f"Valid sequential FG modifications: {valid_fg_mods}")

    return strategy_detected
