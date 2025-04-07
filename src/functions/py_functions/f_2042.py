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
    Detects a convergent synthesis strategy where two complex fragments
    are joined in a late-stage reaction.
    """
    convergent_detected = False

    def calculate_complexity(mol):
        """Calculate molecular complexity based on multiple factors"""
        if mol is None:
            return 0

        atom_count = mol.GetNumAtoms()
        ring_count = mol.GetRingInfo().NumRings()

        # More sophisticated complexity measure
        complexity = atom_count + (ring_count * 3)

        # Add complexity for heteroatoms
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in [1, 6]:  # Not H or C
                complexity += 1

        # Add complexity for functional groups
        complex_fgs = [
            "Ester",
            "Amide",
            "Carboxylic acid",
            "Nitrile",
            "Nitro group",
            "Sulfonamide",
            "Sulfone",
            "Phosphate ester",
        ]
        for fg in complex_fgs:
            if checker.check_fg(fg, Chem.MolToSmiles(mol)):
                complexity += 2

        return complexity

    def is_bond_forming_reaction(rsmi):
        """Check if the reaction forms new carbon-carbon or carbon-heteroatom bonds"""
        # Common bond-forming reaction types
        bond_forming_reactions = [
            "Suzuki",
            "Negishi",
            "Stille",
            "Sonogashira",
            "Heck",
            "Buchwald-Hartwig",
            "N-arylation",
            "Acylation",
            "Mitsunobu",
            "Williamson Ether Synthesis",
            "Ullmann",
            "Goldberg",
            "Michael addition",
            "Aldol",
            "Wittig",
            "Grignard",
            "Diels-Alder",
        ]

        for rxn_type in bond_forming_reactions:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected bond-forming reaction: {rxn_type}")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_detected

        if node["type"] == "reaction" and depth <= 3:  # Expanded to depth 3
            print(f"Examining reaction at depth {depth}")

            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found")
                return

            print(f"Reaction SMILES: {rsmi}")
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least two reactants
            if len(reactants_smiles) >= 2:
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Calculate complexity for each reactant
                    complexities = [calculate_complexity(mol) for mol in reactant_mols]
                    product_complexity = calculate_complexity(product_mol)

                    print(f"Reactant complexities: {complexities}")
                    print(f"Product complexity: {product_complexity}")

                    # Count significant reactants (complexity > 10)
                    significant_reactants = sum(1 for c in complexities if c > 10)

                    # Check for common coupling reactions
                    is_coupling = is_bond_forming_reaction(rsmi)

                    # Check for functional groups involved in coupling
                    coupling_fgs = [
                        "Carboxylic acid",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Boronic acid",
                        "Boronic ester",
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Alkyne",
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Acyl halide",
                        "Ester",
                        "Anhydride",
                    ]

                    fg_count = 0
                    for r_smiles in reactants_smiles:
                        for fg in coupling_fgs:
                            if checker.check_fg(fg, r_smiles):
                                fg_count += 1
                                print(f"Found {fg} in reactant")
                                break

                    # Check for ring structures that might indicate complex fragments
                    ring_structures = [
                        "benzene",
                        "pyridine",
                        "pyrimidine",
                        "pyrazine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "indole",
                        "benzimidazole",
                        "quinoline",
                        "isoquinoline",
                    ]

                    ring_count = 0
                    for r_smiles in reactants_smiles:
                        for ring in ring_structures:
                            if checker.check_ring(ring, r_smiles):
                                ring_count += 1
                                print(f"Found {ring} ring in reactant")
                                break

                    # Criteria for convergent synthesis:
                    # 1. At least 2 significant reactants
                    # 2. Either a known coupling reaction or presence of coupling functional groups
                    # 3. Product should be more complex than individual reactants
                    # 4. Reactants should have some complexity balance

                    if (
                        significant_reactants >= 2
                        and (is_coupling or fg_count >= 2 or ring_count >= 2)
                        and product_complexity > max(complexities)
                    ):

                        # Check relative complexity - reactants should be somewhat balanced
                        sorted_complexities = sorted(complexities, reverse=True)
                        if (
                            len(sorted_complexities) >= 2
                            and sorted_complexities[1] > sorted_complexities[0] * 0.3
                        ):
                            convergent_detected = True
                            print(f"✓ Detected convergent synthesis at depth {depth}")
                        else:
                            print(
                                "Reactant complexities too imbalanced for convergent synthesis"
                            )
                    else:
                        if significant_reactants < 2:
                            print("Not enough significant reactants")
                        if not (is_coupling or fg_count >= 2 or ring_count >= 2):
                            print(
                                "No coupling reaction or sufficient functional groups detected"
                            )
                        if product_complexity <= max(complexities):
                            print("Product not more complex than reactants")

                except Exception as e:
                    print(f"Error in convergent synthesis detection: {e}")

            # Check if this is a hydrolysis or deprotection reaction (not convergent)
            hydrolysis_reactions = [
                "Hydrolysis",
                "Ester saponification",
                "deprotection",
            ]
            is_hydrolysis = any(
                checker.check_reaction(rxn, rsmi) for rxn in hydrolysis_reactions
            )
            if is_hydrolysis:
                print("Detected hydrolysis or deprotection reaction (not convergent)")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Convergent synthesis detected: {convergent_detected}")

    return convergent_detected
