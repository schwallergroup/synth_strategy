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


def main(route, min_count=3):
    """
    This function detects if the synthesis involves multiple C-C bond formations.
    It checks both reaction types and actual structural changes to verify C-C bond formation.
    """
    cc_bond_count = 0

    # List of reactions that typically form C-C bonds
    cc_bond_forming_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with sulfonic esters",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with boronic esters",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Oxidative Heck reaction",
        "Oxidative Heck reaction with vinyl ester",
        "Heck reaction with vinyl ester and amine",
        "Negishi coupling",
        "Stille reaction_vinyl",
        "Stille reaction_aryl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_vinyl OTf",
        "Stille reaction_aryl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "Grignard with CO2 to carboxylic acid",
        "Grignard from aldehyde to alcohol",
        "Grignard from ketone to alcohol",
        "Grignard from nitrile to ketone",
        "Olefination of ketones with Grignard reagents",
        "Olefination of aldehydes with Grignard reagents",
        "Wittig reaction with triphenylphosphorane",
        "Wittig with Phosphonium",
        "Julia Olefination",
        "Aldol condensation",
        "Knoevenagel Condensation",
        "Michael addition",
        "Michael addition methyl",
        "Diels-Alder",
        "Friedel-Crafts alkylation",
        "Friedel-Crafts alkylation with halide",
        "Friedel-Crafts acylation",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl OTf",
        "Sonogashira alkyne_aryl OTf",
        "Sonogashira acetylene_alkenyl halide",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl OTf",
        "Sonogashira alkyne_alkenyl OTf",
        "Sonogashira acetylene_acyl halide",
        "Sonogashira alkyne_acyl halide",
        "Catellani reaction ortho",
        "Catellani reaction para",
        "beta C(sp3) arylation",
        "decarboxylative_coupling",
        "Csp3–Csp2 cross-coupling of alkylarenes to aldehydes",
        "Coupling of benzylboronic esters with aldehydes",
        "Coupling of benzylboronic esters with imines",
        "Pauson-Khand reaction",
        "A3 coupling",
        "A3 coupling to imidazoles",
        "Alkyne-imine cycloaddition",
        "Homologation of aldehydes with formaldehyde",
        "Homologation of diazo compounds to aldehydes using formaldehyde",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Henry Reaction",
        # Additional common C-C bond forming reactions
        "Suzuki",
        "Heck",
        "Stille",
        "Negishi",
        "Wittig",
        "Grignard",
        "Aldol",
        "Michael",
        "Diels-Alder",
        "Friedel-Crafts",
        "Sonogashira",
        "Knoevenagel",
        "decarboxylative coupling",
        "Horner-Wadsworth-Emmons",  # Added this common C-C bond forming reaction
        "Formylation",  # Added this common C-C bond forming reaction
    ]

    def has_new_cc_bond(reactants_smiles, product_smiles, rsmi):
        """Check if a new C-C bond is formed between reactants and product"""
        try:
            # Check for Grignard reactions specifically
            for r in reactants_smiles:
                if "Mg" in r and "[CH" in r:
                    print(f"  Detected potential Grignard reagent: {r}")
                    for other_r in reactants_smiles:
                        if "C=" in other_r and "Mg" not in other_r:
                            print(f"  Detected potential carbonyl compound: {other_r}")
                            print(f"  Likely Grignard C-C bond formation")
                            return True

            # Convert SMILES to molecules
            reactants_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if r.strip()
            ]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if not all(reactants_mols) or not product_mol:
                print(
                    f"Warning: Could not parse some molecules: {reactants_smiles} -> {product_smiles}"
                )
                return False

            # Look for atom-mapped carbons in the product that form new C-C bonds
            # This is more accurate than just counting total bonds
            atom_map_to_idx = {}
            for atom in product_mol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num > 0:
                    atom_map_to_idx[map_num] = atom.GetIdx()

            # Find new C-C bonds in the product that weren't in the reactants
            new_cc_bonds = []
            for bond in product_mol.GetBonds():
                begin_atom = product_mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                end_atom = product_mol.GetAtomWithIdx(bond.GetEndAtomIdx())

                if begin_atom.GetSymbol() == "C" and end_atom.GetSymbol() == "C":
                    begin_map = begin_atom.GetAtomMapNum()
                    end_map = end_atom.GetAtomMapNum()

                    # If both atoms have map numbers, check if they were connected in reactants
                    if begin_map > 0 and end_map > 0:
                        # Try to find these mapped atoms in reactants
                        atoms_connected_in_reactants = False

                        for reactant_mol in reactants_mols:
                            reactant_map_to_idx = {}
                            for atom in reactant_mol.GetAtoms():
                                map_num = atom.GetAtomMapNum()
                                if map_num > 0:
                                    reactant_map_to_idx[map_num] = atom.GetIdx()

                            # Check if both mapped atoms exist in this reactant
                            if (
                                begin_map in reactant_map_to_idx
                                and end_map in reactant_map_to_idx
                            ):
                                begin_idx = reactant_map_to_idx[begin_map]
                                end_idx = reactant_map_to_idx[end_map]

                                # Check if they're bonded in the reactant
                                bond = reactant_mol.GetBondBetweenAtoms(
                                    begin_idx, end_idx
                                )
                                if bond is not None:
                                    atoms_connected_in_reactants = True
                                    break

                        if not atoms_connected_in_reactants:
                            new_cc_bonds.append((begin_map, end_map))

            if new_cc_bonds:
                print(f"  Found {len(new_cc_bonds)} new C-C bonds: {new_cc_bonds}")
                return True

            # Fallback to simpler method if atom mapping doesn't work
            # Count C-C bonds in reactants and product
            def count_cc_bonds(mol):
                count = 0
                for bond in mol.GetBonds():
                    begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
                    end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
                    if begin_atom.GetSymbol() == "C" and end_atom.GetSymbol() == "C":
                        count += 1
                return count

            reactants_cc_bonds = sum(count_cc_bonds(mol) for mol in reactants_mols)
            product_cc_bonds = count_cc_bonds(product_mol)

            # If product has more C-C bonds than reactants combined, a new C-C bond was formed
            if product_cc_bonds > reactants_cc_bonds:
                print(
                    f"  Fallback method: Product has {product_cc_bonds} C-C bonds, reactants have {reactants_cc_bonds}"
                )
                return True

            return False
        except Exception as e:
            print(f"Error in has_new_cc_bond: {e}")
            return False

    def is_formylation_reaction(rsmi):
        """Check if the reaction is a formylation (adding a formyl group)"""
        try:
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains formaldehyde or a formylating agent
            formylating_agents = ["C=O", "CH2O", "CHO", "N(C)C=O"]
            has_formylating_agent = any(
                any(agent in r for agent in formylating_agents) for r in reactants
            )

            # Check if the product has a new formyl group
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                for atom in product_mol.GetAtoms():
                    if atom.GetSymbol() == "C":
                        neighbors = atom.GetNeighbors()
                        if any(
                            n.GetSymbol() == "O"
                            and n.GetTotalNumHs() == 0
                            and n.GetDegree() == 1
                            for n in neighbors
                        ):
                            if has_formylating_agent:
                                print(f"  Detected formylation reaction")
                                return True

            return False
        except Exception as e:
            print(f"Error in is_formylation_reaction: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal cc_bond_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            try:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants_part = parts[0]
                    reagents_part = parts[1]
                    product_part = parts[2]

                    # Split reactants by "." but handle empty strings
                    reactants = [r for r in reactants_part.split(".") if r.strip()]
                    product = product_part

                    # Method 1: Check if this reaction is a known C-C bond forming reaction
                    reaction_match = False
                    for reaction_type in cc_bond_forming_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"  Matched known C-C bond forming reaction: {reaction_type}"
                            )
                            reaction_match = True
                            break

                    # Check for Horner-Wadsworth-Emmons reaction
                    if not reaction_match and "P(=O)" in rsmi and "C=" in product:
                        print(f"  Detected potential Horner-Wadsworth-Emmons reaction")
                        reaction_match = True

                    # Check for formylation
                    if not reaction_match and is_formylation_reaction(rsmi):
                        reaction_match = True

                    # Method 2: Check for actual C-C bond formation by comparing structures
                    structure_match = has_new_cc_bond(reactants, product, rsmi)
                    if structure_match:
                        print(f"  Structural analysis confirms new C-C bond formation")

                    # Count if either method confirms C-C bond formation
                    if reaction_match or structure_match:
                        cc_bond_count += 1
                        print(
                            f"  C-C bond formation detected! Current count: {cc_bond_count}"
                        )

            except Exception as e:
                print(f"Error processing reaction {rsmi}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    print("Starting analysis for multiple C-C bond formations...")
    dfs_traverse(route)

    result = cc_bond_count >= min_count
    print(f"Multiple C-C bond formations ({cc_bond_count} >= {min_count}): {result}")
    return result
