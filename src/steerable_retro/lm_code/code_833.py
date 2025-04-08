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
    This function detects if the synthetic route involves a ring opening strategy.
    It looks for reactions where a ring is opened during the transformation.
    """
    ring_opening_detected = False

    def dfs_traverse(node):
        nonlocal ring_opening_detected

        if ring_opening_detected:
            return  # Early return if already detected

        if node["type"] == "reaction":
            # Check if RingBreaker flag is set in metadata
            if node.get("metadata", {}).get("RingBreaker", False):
                print(f"Ring opening detected via RingBreaker flag in metadata")
                ring_opening_detected = True
                return

            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for known ring opening reactions
                ring_opening_reactions = [
                    "Acetal hydrolysis to diol",
                    "Acetal hydrolysis to aldehyde",
                    "Ketal hydrolysis to ketone",
                    "Retro-Diels-Alder from oxazole",
                    "Ring opening of epoxide with amine",
                ]

                for rxn_type in ring_opening_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Ring opening detected: {rxn_type} in {rsmi}")
                        ring_opening_detected = True
                        return

                # If no specific reaction type matched, check ring counts
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Create molecules with error handling
                reactant_mols = []
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        reactant_mols.append(mol)

                product_mol = Chem.MolFromSmiles(product_smiles)
                if not product_mol or not reactant_mols:
                    return

                # Count rings in each reactant and the product
                reactants_ring_info = [(mol, mol.GetRingInfo().NumRings()) for mol in reactant_mols]
                product_ring_count = product_mol.GetRingInfo().NumRings()

                # Check for common ring structures that might be opened
                ring_types = [
                    "oxirane",
                    "aziridine",
                    "cyclopropane",
                    "oxetane",
                    "azetidine",
                    "cyclobutane",
                    "tetrahydrofuran",
                    "pyrrolidine",
                    "cyclopentane",
                    "tetrahydropyran",
                    "piperidine",
                    "cyclohexane",
                ]

                # Check if any reactant has more rings than the product
                for mol, ring_count in reactants_ring_info:
                    if ring_count > 0 and product_ring_count < ring_count:
                        # Check if any specific ring type was present in reactant but not in product
                        for ring_type in ring_types:
                            if checker.check_ring(
                                ring_type, Chem.MolToSmiles(mol)
                            ) and not checker.check_ring(ring_type, product_smiles):
                                print(
                                    f"Ring opening detected: {ring_type} ring in reactant but not in product: {rsmi}"
                                )
                                ring_opening_detected = True
                                return

                # Check total ring count as a fallback
                total_reactant_rings = sum(count for _, count in reactants_ring_info)
                if product_ring_count < total_reactant_rings:
                    # This is a potential ring opening, but we need to verify it's not just
                    # a reactant with rings that wasn't incorporated

                    # If there's only one reactant, it's definitely ring opening
                    if len(reactant_mols) == 1:
                        print(f"Ring opening detected (single reactant): {rsmi}")
                        ring_opening_detected = True
                    else:
                        # For multiple reactants, check if the product contains atoms from all reactants
                        # This is a simplified check - in a real scenario, you'd want to use atom mapping
                        # to track exactly which atoms from which reactant end up in the product
                        print(f"Potential ring opening detected (multiple reactants): {rsmi}")
                        ring_opening_detected = True

            except Exception as e:
                print(f"Error in ring opening detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return ring_opening_detected
