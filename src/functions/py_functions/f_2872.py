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
    Detects if the synthesis follows a convergent approach with late-stage fragment coupling.
    Specifically looks for a reaction in the first half of the synthesis that combines two complex fragments.
    """
    fragment_coupling_detected = False
    max_depth = 0
    target_mol = None

    def is_coupling_reaction(rsmi):
        """Check if the reaction is a known coupling reaction type"""
        coupling_reactions = [
            # C-C coupling reactions
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic acids OTf",
            "Suzuki coupling with boronic esters",
            "Suzuki coupling with boronic esters OTf",
            "Suzuki coupling with sulfonic esters",
            "Negishi coupling",
            "Stille reaction_aryl",
            "Stille reaction_vinyl",
            "Stille reaction_benzyl",
            "Stille reaction_allyl",
            "Stille reaction_aryl OTf",
            "Stille reaction_vinyl OTf",
            "Stille reaction_benzyl OTf",
            "Stille reaction_allyl OTf",
            "Stille reaction_other",
            "Stille reaction_other OTf",
            "Heck_terminal_vinyl",
            "Heck_non-terminal_vinyl",
            "Oxidative Heck reaction",
            "Oxidative Heck reaction with vinyl ester",
            "Heck reaction with vinyl ester and amine",
            "Sonogashira alkyne_aryl halide",
            "Sonogashira acetylene_aryl halide",
            "Sonogashira alkyne_aryl OTf",
            "Sonogashira acetylene_aryl OTf",
            "Sonogashira alkyne_alkenyl halide",
            "Sonogashira acetylene_alkenyl halide",
            "Sonogashira alkyne_alkenyl OTf",
            "Sonogashira acetylene_alkenyl OTf",
            "Sonogashira alkyne_acyl halide",
            "Sonogashira acetylene_acyl halide",
            "Kumada cross-coupling",
            "Hiyama-Denmark Coupling",
            "decarboxylative_coupling",
            "Catellani reaction ortho",
            "Catellani reaction para",
            "Aryllithium cross-coupling",
            # C-N coupling reactions
            "Buchwald-Hartwig",
            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
            "Goldberg coupling",
            "Goldberg coupling aryl amine-aryl chloride",
            "Goldberg coupling aryl amide-aryl chloride",
            "Ullmann-Goldberg Substitution amine",
            "N-arylation_heterocycles",
            # C-O coupling reactions
            "Ullmann-Goldberg Substitution aryl alcohol",
            "Ullmann-Goldberg Substitution thiol",
            "Ullmann condensation",
            "Chan-Lam alcohol",
            "Chan-Lam amine",
            "Chan-Lam etherification",
        ]

        for rxn_type in coupling_reactions:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected coupling reaction: {rxn_type}")
                return True
        return False

    def is_complex_fragment(smiles):
        """Determine if a fragment is complex based on multiple criteria"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check size and rings
        atom_count = mol.GetNumAtoms()
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()

        # Count heteroatoms (non-carbon, non-hydrogen)
        heteroatom_count = sum(
            1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6]
        )

        # Check for functional groups that indicate complexity
        complex_fgs = [
            "Aromatic halide",
            "Boronic acid",
            "Boronic ester",
            "Ester",
            "Amide",
            "Nitrile",
            "Ether",
            "Alcohol",
            "Carboxylic acid",
            "Phenol",
            "Aniline",
            "Ketone",
            "Aldehyde",
            "Sulfonamide",
            "Sulfone",
            "Phosphate ester",
            "Triflate",
            "Tosylate",
            "Mesylate",
            "Azide",
            "Nitro group",
            "Primary amine",
            "Secondary amine",
            "Tertiary amine",
        ]

        has_complex_fg = any(checker.check_fg(fg, smiles) for fg in complex_fgs)

        # Check for complex ring systems
        complex_rings = [
            "indole",
            "quinoline",
            "isoquinoline",
            "benzothiophene",
            "benzofuran",
            "benzimidazole",
            "naphthalene",
            "pyridine",
            "pyrimidine",
            "pyrazine",
            "piperidine",
            "morpholine",
            "furan",
            "thiophene",
            "pyrrole",
            "imidazole",
            "oxazole",
            "thiazole",
            "triazole",
            "tetrazole",
            "benzene",
        ]

        has_complex_ring = any(
            checker.check_ring(ring, smiles) for ring in complex_rings
        )

        # A fragment is complex if it meets multiple criteria
        complexity_score = 0
        if atom_count >= 7:
            complexity_score += 1
        if ring_count >= 1:
            complexity_score += 1
        if heteroatom_count >= 2:
            complexity_score += 1
        if has_complex_fg:
            complexity_score += 1
        if has_complex_ring:
            complexity_score += 1

        return complexity_score >= 2

    def get_fragment_complexity_score(smiles):
        """Get a numerical complexity score for a fragment"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0

        score = 0
        # Size contribution
        atom_count = mol.GetNumAtoms()
        score += min(5, atom_count // 3)  # Cap size contribution, more sensitive

        # Ring contribution
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        score += min(4, ring_count * 2)  # Increased weight for rings

        # Heteroatom contribution
        heteroatom_count = sum(
            1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6]
        )
        score += min(3, heteroatom_count)  # More weight to heteroatoms

        # Functional group contribution
        complex_fgs = [
            "Aromatic halide",
            "Boronic acid",
            "Boronic ester",
            "Ester",
            "Amide",
            "Nitrile",
            "Ether",
            "Alcohol",
            "Carboxylic acid",
            "Phenol",
            "Aniline",
            "Ketone",
            "Aldehyde",
            "Sulfonamide",
            "Sulfone",
        ]

        fg_count = sum(1 for fg in complex_fgs if checker.check_fg(fg, smiles))
        score += min(3, fg_count)

        return score

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_detected, max_depth, target_mol

        # Store the target molecule (root of the tree)
        if depth == 0 and node["type"] == "mol":
            target_mol = Chem.MolFromSmiles(node["smiles"])

        # Check for late-stage coupling reactions
        # Consider "late-stage" as first 2/3 of the synthesis depth
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we're combining multiple fragments
            if len(reactants) >= 2:
                # Count complex fragments
                complex_fragments = [r for r in reactants if is_complex_fragment(r)]

                # For late-stage reactions (closer to target)
                if depth <= max(3, max_depth * 2 // 3):
                    # Check if it's a coupling reaction or if we have exactly two complex fragments
                    if len(complex_fragments) >= 2 and (
                        is_coupling_reaction(rsmi) or len(reactants) == 2
                    ):
                        print(f"Late-stage fragment coupling detected at depth {depth}")
                        print(f"Reactants: {reactants}")

                        # Check relative sizes of fragments to ensure balanced convergence
                        complexity_scores = [
                            get_fragment_complexity_score(r) for r in complex_fragments
                        ]
                        min_score = min(complexity_scores)
                        max_score = max(complexity_scores)

                        # Ensure fragments are of comparable complexity (smaller one at least 1/4 of larger)
                        if (
                            min_score > 0
                            and max_score > 0
                            and min_score >= max_score / 4
                        ):
                            print(f"Fragment complexity scores: {complexity_scores}")
                            print(
                                f"Balanced fragment coupling detected (min/max ratio: {min_score/max_score:.2f})"
                            )
                            fragment_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First pass to determine max_depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")

    # Second pass to detect coupling
    dfs_traverse(route)

    # Only consider it convergent if we have late-stage coupling and sufficient depth
    is_convergent = fragment_coupling_detected and max_depth >= 3
    print(f"Convergent synthesis with late-stage coupling: {is_convergent}")
    return is_convergent
