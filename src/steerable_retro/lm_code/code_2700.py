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
    Detects if the route involves formation of a heterocycle containing a trifluoromethyl group.
    """
    trifluoromethyl_heterocycle_formed = False

    # List of heterocycles to check
    heterocycle_types = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    def is_cf3_attached_to_heterocycle(mol_smiles):
        """Check if a CF3 group is attached to a heterocycle in the molecule."""
        mol = Chem.MolFromSmiles(mol_smiles)
        if mol is None:
            return False

        # First check if molecule has a trifluoro group
        if not checker.check_fg("Trifluoro group", mol_smiles):
            return False

        # Then check for heterocycles
        for heterocycle in heterocycle_types:
            if checker.check_ring(heterocycle, mol_smiles):
                # Get atom indices for both the CF3 group and the heterocycle
                try:
                    cf3_indices = checker.get_fg_atom_indices("Trifluoro group", mol_smiles)
                    ring_indices = checker.get_ring_atom_indices(heterocycle, mol_smiles)

                    if not cf3_indices or not ring_indices:
                        continue

                    # Process CF3 indices to get carbon atoms
                    cf3_carbons = []
                    for cf3_group in cf3_indices:
                        # Handle different possible return structures
                        if isinstance(cf3_group, tuple) and len(cf3_group) > 0:
                            if isinstance(cf3_group[0], tuple) and len(cf3_group[0]) > 0:
                                cf3_carbons.append(cf3_group[0][0])  # Nested tuple
                            else:
                                cf3_carbons.append(cf3_group[0])  # Simple tuple
                        elif isinstance(cf3_group, int):
                            cf3_carbons.append(cf3_group)  # Direct integer

                    # Process ring indices to get all ring atoms
                    all_ring_atoms = set()
                    for ring_instance in ring_indices:
                        # Handle different possible return structures
                        if isinstance(ring_instance, tuple):
                            if all(isinstance(x, tuple) for x in ring_instance if x):
                                # Nested tuples
                                for atom_tuple in ring_instance:
                                    all_ring_atoms.update(atom_tuple)
                            else:
                                # Simple tuple
                                all_ring_atoms.update(ring_instance)
                        elif isinstance(ring_instance, int):
                            all_ring_atoms.add(ring_instance)  # Direct integer

                    # Check if any CF3 carbon is directly attached to any ring atom
                    for cf3_carbon in cf3_carbons:
                        for ring_atom in all_ring_atoms:
                            bond = mol.GetBondBetweenAtoms(cf3_carbon, ring_atom)
                            if bond is not None:
                                print(
                                    f"Found CF3-heterocycle attachment: CF3 carbon {cf3_carbon} is bonded to {heterocycle} atom {ring_atom}"
                                )
                                return True

                except Exception as e:
                    print(f"Error checking CF3-heterocycle attachment for {heterocycle}: {e}")
                    print(f"CF3 indices structure: {type(cf3_indices)}, content: {cf3_indices}")
                    if cf3_indices:
                        print(f"First CF3 group type: {type(cf3_indices[0])}")
                    print(f"Ring indices structure: {type(ring_indices)}, content: {ring_indices}")
                    if ring_indices:
                        print(f"First ring instance type: {type(ring_indices[0])}")
                    continue

        return False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_heterocycle_formed

        if trifluoromethyl_heterocycle_formed:
            return  # Early return if already found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if product has a CF3-heterocycle
                product_has_cf3_heterocycle = is_cf3_attached_to_heterocycle(product_smiles)

                if product_has_cf3_heterocycle:
                    print(f"Product has CF3-heterocycle: {product_smiles}")
                    # Check if any reactant already has the CF3-heterocycle
                    reactants_with_cf3 = []
                    for r in reactants_smiles:
                        if is_cf3_attached_to_heterocycle(r):
                            reactants_with_cf3.append(r)

                    if not reactants_with_cf3:
                        print(f"Detected trifluoromethyl heterocycle formation at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        trifluoromethyl_heterocycle_formed = True
                        return
                    else:
                        print(f"CF3-heterocycle already present in reactants: {reactants_with_cf3}")
            except Exception as e:
                print(f"Error processing reaction node at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return trifluoromethyl_heterocycle_formed
