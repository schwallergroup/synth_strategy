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
    This function detects a strategy involving sequential functionalization of an aromatic ring.
    """
    # Track functionalization steps on aromatic rings
    # Store (reaction_type, depth, product_smiles, ring_atoms, modified_atom_indices) tuples
    aromatic_modifications = []

    # List of functional groups to check for aromatic functionalization
    aromatic_fg_types = [
        "Nitro group",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aromatic halide",
        "Phenol",
        "Aromatic alcohol",
        "Triflate",
        "Tosylate",
        "Mesylate",
        "Sulfonamide",
        "Sulfonic acid",
        "Sulfonate",
        "Nitrile",
        "Carboxylic acid",
        "Ester",
        "Aldehyde",
        "Ketone",
        "Amide",
    ]

    # List of reaction types that typically involve aromatic functionalization
    aromatic_reaction_types = [
        "Aromatic nitration with HNO3",
        "Aromatic nitration with NO3 salt",
        "Aromatic nitration with NO2 salt",
        "Aromatic nitration with alkyl NO2",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "Friedel-Crafts alkylation",
        "Friedel-Crafts acylation",
        "Suzuki coupling with boronic acids",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Aromatic hydroxylation",
        "Aromatic sulfonyl chlorination",
        "Minisci (para)",
        "Minisci (ortho)",
    ]

    # List of aromatic rings to check
    aromatic_rings = [
        "benzene",
        "pyridine",
        "furan",
        "thiophene",
        "pyrrole",
        "naphthalene",
        "anthracene",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is an aromatic functionalization reaction
            is_aromatic_functionalization = False
            reaction_type = None
            ring_atoms = None
            modified_atoms = None

            # First check if it matches any known aromatic functionalization reaction type
            for rxn_type in aromatic_reaction_types:
                if checker.check_reaction(rxn_type, rsmi):
                    is_aromatic_functionalization = True
                    reaction_type = rxn_type
                    print(f"Found {rxn_type} at depth {depth}")

                    # Find which aromatic ring was modified
                    for ring_type in aromatic_rings:
                        if checker.check_ring(ring_type, product_smiles):
                            ring_atoms = checker.get_ring_atom_indices(ring_type, product_smiles)
                            if ring_atoms:
                                break
                    break

            # If no specific reaction type matched, check if it involves adding a functional group to an aromatic ring
            if not is_aromatic_functionalization:
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None
                if product_mol:
                    # Check if product has an aromatic ring
                    has_aromatic_ring = False
                    for ring_type in aromatic_rings:
                        if checker.check_ring(ring_type, product_smiles):
                            has_aromatic_ring = True
                            ring_atoms = checker.get_ring_atom_indices(ring_type, product_smiles)
                            break

                    if has_aromatic_ring:
                        # Check for functional groups on the aromatic ring
                        for fg_type in aromatic_fg_types:
                            if checker.check_fg(fg_type, product_smiles):
                                # Get the atoms involved in this functional group
                                fg_atoms = checker.get_fg_atom_indices(fg_type, product_smiles)

                                # Check if this FG was added in this reaction
                                fg_in_product = True
                                fg_in_any_reactant = False

                                for reactant in reactants_smiles:
                                    if checker.check_fg(fg_type, reactant):
                                        fg_in_any_reactant = True
                                        break

                                if fg_in_product and not fg_in_any_reactant:
                                    is_aromatic_functionalization = True
                                    reaction_type = f"Addition of {fg_type}"
                                    modified_atoms = fg_atoms
                                    print(f"Found addition of {fg_type} at depth {depth}")
                                    break

            if is_aromatic_functionalization and reaction_type and ring_atoms:
                aromatic_modifications.append(
                    (reaction_type, depth, product_smiles, ring_atoms, modified_atoms)
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort modifications by depth (to get chronological order in synthesis direction)
    aromatic_modifications.sort(key=lambda x: x[1], reverse=True)

    # Check if we have at least 2 modifications
    if len(aromatic_modifications) >= 2:
        print("Detected potential sequential aromatic functionalization strategy")
        for mod_type, depth, _, _, _ in aromatic_modifications:
            print(f"  - {mod_type} at depth {depth}")

        # Check if modifications are sequential (on the same ring)
        # Allow for reasonable gaps between modifications (up to 5 steps)
        sequential = True
        for i in range(len(aromatic_modifications) - 1):
            current_mod = aromatic_modifications[i]
            next_mod = aromatic_modifications[i + 1]

            # Check if the modifications are reasonably close in the synthesis sequence
            if abs(current_mod[1] - next_mod[1]) > 5:  # Allow reasonable gaps between modifications
                sequential = False
                print(
                    f"Too large gap between modifications at depths {current_mod[1]} and {next_mod[1]}"
                )
                break

            # Check if modifications are on the same ring
            # This is a simplified check - in a real implementation, you would need to
            # track atom mappings between reactions to ensure it's the same ring
            current_ring = current_mod[3]
            next_ring = next_mod[3]

            if not current_ring or not next_ring:
                continue  # Skip if we couldn't identify ring atoms

            # If we have at least one ring identified in both steps, consider it sequential
            # In a more sophisticated implementation, you would check atom mappings
            # to ensure it's the same ring across reactions

        return sequential

    return False
