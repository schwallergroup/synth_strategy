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
    Detects a synthetic strategy involving nitration of aromatic rings.
    """
    nitration_count = 0
    aromatic_nitro_molecules = set()  # Track molecules with aromatic nitro groups

    def dfs_traverse(node, depth=0):
        nonlocal nitration_count, aromatic_nitro_molecules

        if node["type"] == "reaction":
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a known aromatic nitration reaction type
            nitration_reaction_types = [
                "Aromatic nitration with HNO3",
                "Aromatic nitration with NO3 salt",
                "Aromatic nitration with NO2 salt",
                "Aromatic nitration with alkyl NO2",
            ]

            is_nitration = False
            for rxn_type in nitration_reaction_types:
                if checker.check_reaction(rxn_type, rsmi):
                    is_nitration = True
                    print(f"Found {rxn_type} reaction at depth {depth}: {rsmi}")
                    break

            # If not identified by reaction type, check for nitro group addition to aromatic ring
            if not is_nitration:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has a nitro group
                if checker.check_fg("Nitro group", product_smiles):
                    # Check if any reactant has an aromatic ring
                    has_aromatic_reactant = False

                    for reactant in reactants_smiles:
                        # Check for aromatic rings in reactants using RDKit
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            for atom in reactant_mol.GetAtoms():
                                if atom.GetIsAromatic():
                                    has_aromatic_reactant = True
                                    break
                            if has_aromatic_reactant:
                                break

                    # Count nitro groups in reactants
                    total_reactant_nitro = 0
                    for reactant in reactants_smiles:
                        reactant_nitro_indices = checker.get_fg_atom_indices(
                            "Nitro group", reactant
                        )
                        if reactant_nitro_indices:
                            total_reactant_nitro += len(reactant_nitro_indices)

                    # Count nitro groups in product
                    product_nitro_indices = checker.get_fg_atom_indices(
                        "Nitro group", product_smiles
                    )
                    product_nitro_count = len(product_nitro_indices) if product_nitro_indices else 0

                    # Check if nitro groups are attached to aromatic atoms in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_nitro_indices:
                        aromatic_nitro_count = 0
                        for nitro_group in product_nitro_indices:
                            # The first atom in the nitro group pattern is the one attached to the ring
                            nitro_attachment_idx = nitro_group[0]
                            if nitro_attachment_idx < product_mol.GetNumAtoms():
                                attachment_atom = product_mol.GetAtomWithIdx(nitro_attachment_idx)
                                for neighbor in attachment_atom.GetNeighbors():
                                    if neighbor.GetIsAromatic():
                                        aromatic_nitro_count += 1
                                        break

                        # If product has more aromatic nitro groups than reactants and there's an aromatic reactant
                        if aromatic_nitro_count > 0 and has_aromatic_reactant:
                            is_nitration = True
                            print(f"Found nitration by nitro group attachment at depth {depth}")

            if is_nitration:
                nitration_count += 1

        elif node["type"] == "mol":
            # Check if this molecule has a nitro group attached to an aromatic ring
            mol_smiles = node["smiles"]

            if checker.check_fg("Nitro group", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    nitro_indices = checker.get_fg_atom_indices("Nitro group", mol_smiles)
                    if nitro_indices:
                        for nitro_group in nitro_indices:
                            nitro_attachment_idx = nitro_group[0]
                            if nitro_attachment_idx < mol.GetNumAtoms():
                                attachment_atom = mol.GetAtomWithIdx(nitro_attachment_idx)
                                for neighbor in attachment_atom.GetNeighbors():
                                    if neighbor.GetIsAromatic():
                                        aromatic_nitro_molecules.add(mol_smiles)
                                        print(
                                            f"Found molecule with aromatic nitro group at depth {depth}: {mol_smiles}"
                                        )
                                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if at least one nitration is detected or if molecules with aromatic nitro groups exist
    is_strategy_present = nitration_count > 0 or len(aromatic_nitro_molecules) > 0

    print(f"Nitration count: {nitration_count}")
    print(f"Molecules with aromatic nitro groups: {len(aromatic_nitro_molecules)}")
    print(f"Aromatic nitration strategy detected: {is_strategy_present}")

    return is_strategy_present
