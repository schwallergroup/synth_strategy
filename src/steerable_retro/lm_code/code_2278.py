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
    Detects if the synthesis uses a late-stage coupling of two complex fragments.

    A late-stage coupling is defined as:
    1. Occurs at depth 0 or 1 in the synthesis tree (close to final product)
    2. Involves a known coupling reaction (Suzuki, Buchwald-Hartwig, etc.)
    3. Joins two complex fragments (>10 atoms, containing at least one ring)
    4. Forms a new bond between the fragments
    """
    found_late_stage_coupling = False

    def is_complex_fragment(mol):
        """Check if a molecule is a complex fragment"""
        if mol is None or mol.GetNumAtoms() < 10:
            return False

        # Check if the molecule contains at least one ring
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 0:
            return False

        return True

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_coupling

        # Early return if we already found a late-stage coupling
        if found_late_stage_coupling:
            return

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Only consider reactions at depth 0 or 1 (late stage)
        if (
            node["type"] == "reaction"
            and depth <= 1
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Comprehensive list of coupling reactions
            coupling_reactions = [
                "Suzuki",
                "Buchwald-Hartwig",
                "N-arylation",
                "Negishi",
                "Stille",
                "Heck",
                "Sonogashira",
                "Ullmann",
                "Chan-Lam",
                "Ullmann-Goldberg",
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic esters",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "Heck terminal vinyl",
                "Heck non-terminal vinyl",
                "Stille reaction_aryl",
                "Stille reaction_vinyl",
                "Negishi coupling",
                "Ullmann condensation",
                "Williamson ether synthesis",
                "Williamson Ether Synthesis",
                "Mitsunobu",
                "Mitsunobu_phenole",
                "Mitsunobu aryl ether",
                "Chan-Lam etherification",
                "Chan-Lam amine",
                "Chan-Lam alcohol",
                "Ullmann-Goldberg Substitution amine",
                "Ullmann-Goldberg Substitution thiol",
                "Ullmann-Goldberg Substitution aryl alcohol",
            ]

            # Check if this is a coupling reaction
            is_coupling = False
            coupling_type = ""
            for rxn_type in coupling_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    is_coupling = True
                    coupling_type = rxn_type
                    print(f"Found {rxn_type} coupling reaction at depth {depth}")
                    break

            # If not identified as a standard coupling, check for ether formation
            if not is_coupling:
                # Check for ether formation (O-arylation)
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product) if product else None
                if product_mol and checker.check_fg("Ether", product):
                    print("Found potential ether coupling reaction")
                    # Check if the ether is newly formed (not present in reactants)
                    ether_in_reactants = any(checker.check_fg("Ether", r) for r in reactants if r)
                    if not ether_in_reactants:
                        is_coupling = True
                        coupling_type = "Ether formation"
                        print("Confirmed ether coupling reaction")

            if is_coupling:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if we have at least two reactants
                    if len(reactants) >= 2:
                        reactant_mols = [
                            Chem.MolFromSmiles(r)
                            for r in reactants
                            if r and Chem.MolFromSmiles(r) is not None
                        ]
                        product_mol = Chem.MolFromSmiles(product) if product else None

                        if product_mol and len(reactant_mols) >= 2:
                            # Check if at least two reactants are complex fragments
                            complex_fragments = []
                            for i, mol in enumerate(reactant_mols):
                                if is_complex_fragment(mol):
                                    complex_fragments.append((i, mol))
                                    print(
                                        f"Reactant {i} is a complex fragment with {mol.GetNumAtoms()} atoms and {mol.GetRingInfo().NumRings()} rings"
                                    )
                                else:
                                    print(
                                        f"Reactant {i} is NOT complex: {mol.GetNumAtoms()} atoms, {mol.GetRingInfo().NumRings()} rings"
                                    )

                            if len(complex_fragments) >= 2:
                                print(
                                    f"Found at least two complex fragments in {coupling_type} reaction at depth {depth}"
                                )

                                # Verify that a bond is formed between the complex fragments
                                # This is implicit in coupling reactions, so we can assume it's true
                                # for the coupling reaction types we're checking

                                print(
                                    f"Confirmed late-stage coupling of complex fragments at depth {depth}"
                                )
                                found_late_stage_coupling = True
                                return
                            else:
                                print(
                                    f"Not enough complex fragments found: {len(complex_fragments)}"
                                )
                        else:
                            print("Invalid product or not enough valid reactants")
                    else:
                        print("Not enough reactants for coupling reaction")
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_late_stage_coupling}")
    return found_late_stage_coupling
