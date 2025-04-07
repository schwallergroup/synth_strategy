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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects convergent synthesis by identifying steps where
    two complex fragments are combined.
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count complex reactants (using refined complexity criteria)
                complex_reactants = 0
                reactant_mols = []
                for r in reactants:
                    mol = Chem.MolFromSmiles(r) if r else None
                    if mol:
                        reactant_mols.append(mol)
                        # Consider both atom count and ring structures
                        atom_count = mol.GetNumAtoms()
                        ring_atom_count = sum(1 for atom in mol.GetAtoms() if atom.IsInRing())

                        if atom_count > 8 or ring_atom_count > 5:
                            complex_reactants += 1

                # If we have at least 2 complex reactants, it's potentially convergent
                if complex_reactants >= 2:
                    # Check if this is a coupling reaction (common in convergent synthesis)
                    is_coupling = False

                    # List of coupling reactions from the provided reaction types
                    coupling_reactions = [
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
                        "Heck terminal vinyl",
                        "Heck_terminal_vinyl",
                        "Heck_non-terminal_vinyl",
                        "Sonogashira alkyne_aryl halide",
                        "Sonogashira acetylene_aryl halide",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Ullmann-Goldberg Substitution amine",
                        "Ullmann-Goldberg Substitution thiol",
                        "Ullmann-Goldberg Substitution aryl alcohol",
                        "Ullmann condensation",
                        "Hiyama-Denmark Coupling",
                        "Kumada cross-coupling",
                        "Aryllithium cross-coupling",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_coupling = True
                            print(f"Coupling reaction detected: {rxn_type}")
                            break

                    # Also consider other reactions that join fragments
                    fragment_joining_reactions = [
                        "Williamson Ether Synthesis",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Esterification of Carboxylic Acids",
                        "Ugi reaction",
                        "Diels-Alder",
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Mitsunobu esterification",
                        "Mitsunobu aryl ether",
                    ]

                    if not is_coupling:
                        for rxn_type in fragment_joining_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                is_coupling = True
                                print(f"Fragment-joining reaction detected: {rxn_type}")
                                break

                    # Prioritize late-stage convergent steps (lower depth)
                    significance_factor = max(1, 10 - depth)  # Higher for late-stage reactions

                    # Calculate complexity score based on reactant structures
                    complexity_score = 0
                    for mol in reactant_mols:
                        if mol:
                            atom_count = mol.GetNumAtoms()
                            ring_count = len(Chem.GetSSSR(mol))
                            complexity_score += atom_count + (ring_count * 2)

                    # If it's a coupling reaction or we have very complex fragments
                    if is_coupling or complexity_score > 40:
                        print(f"Convergent synthesis detected at depth {depth}: {rsmi}")
                        print(
                            f"Complex reactants: {complex_reactants}, Complexity score: {complexity_score}"
                        )
                        convergent_synthesis_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if not convergent_synthesis_detected:
        print("No convergent synthesis strategy detected in this route")

    return convergent_synthesis_detected
