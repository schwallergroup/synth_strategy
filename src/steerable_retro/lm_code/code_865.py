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
    This function detects preservation of key heterocyclic systems (triazole and pyridone)
    throughout the synthetic route.
    """
    # Track heterocycles through the synthesis
    heterocycle_tracking = {
        "triazole": {"preserved": True, "present_in_final": False},
        "pyridone": {"preserved": True, "present_in_final": False},
    }

    # First pass: identify the final product (root node)
    final_product_smiles = route["smiles"] if route["type"] == "mol" else None
    if not final_product_smiles:
        print("Could not identify final product")
        return False

    print(f"Final product identified: {final_product_smiles}")

    # Check if triazole is present in final product
    if checker.check_ring("triazole", final_product_smiles):
        print(f"Found triazole in final product: {final_product_smiles}")
        heterocycle_tracking["triazole"]["present_in_final"] = True

    # Check if pyridone is present in final product
    pyridone_in_final = False
    if checker.check_ring("pyridine", final_product_smiles) and checker.check_fg(
        "Ketone", final_product_smiles
    ):
        # Get the atom indices for both
        pyridine_indices = checker.get_ring_atom_indices("pyridine", final_product_smiles)
        ketone_indices = checker.get_fg_atom_indices("Ketone", final_product_smiles)

        # Check if any ketone is attached to the pyridine ring
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol:
            for ring_atoms in pyridine_indices:
                for ketone_atoms in ketone_indices:
                    # Check if the carbonyl carbon is connected to the pyridine ring
                    carbonyl_carbon = ketone_atoms[0][0]
                    for bond in mol.GetBonds():
                        begin_idx = bond.GetBeginAtomIdx()
                        end_idx = bond.GetEndAtomIdx()
                        if (begin_idx == carbonyl_carbon and end_idx in ring_atoms[0]) or (
                            end_idx == carbonyl_carbon and begin_idx in ring_atoms[0]
                        ):
                            print(f"Found pyridone in final product: {final_product_smiles}")
                            heterocycle_tracking["pyridone"]["present_in_final"] = True
                            pyridone_in_final = True
                            break
                    if pyridone_in_final:
                        break
                if pyridone_in_final:
                    break

    # Second pass: track heterocycles through reactions
    def track_heterocycles(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if triazole is preserved in this reaction (retrosynthetic perspective)
            if heterocycle_tracking["triazole"]["preserved"]:
                reactant_has_triazole = any(checker.check_ring("triazole", r) for r in reactants)
                product_has_triazole = checker.check_ring("triazole", product)

                # In retrosynthesis, if product has triazole but reactants don't, it means
                # triazole was created in the forward direction (not preserved)
                if product_has_triazole and not reactant_has_triazole:
                    print(f"Triazole not preserved in reaction (created): {rsmi}")
                    heterocycle_tracking["triazole"]["preserved"] = False

            # Check if pyridone is preserved in this reaction
            if heterocycle_tracking["pyridone"]["preserved"]:
                # Check for pyridone in reactants
                reactant_has_pyridone = False
                for r in reactants:
                    if checker.check_ring("pyridine", r) and checker.check_fg("Ketone", r):
                        # Verify it's actually a pyridone by checking connections
                        pyridine_indices = checker.get_ring_atom_indices("pyridine", r)
                        ketone_indices = checker.get_fg_atom_indices("Ketone", r)
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            for ring_atoms in pyridine_indices:
                                for ketone_atoms in ketone_indices:
                                    carbonyl_carbon = ketone_atoms[0][0]
                                    for bond in mol.GetBonds():
                                        begin_idx = bond.GetBeginAtomIdx()
                                        end_idx = bond.GetEndAtomIdx()
                                        if (
                                            begin_idx == carbonyl_carbon
                                            and end_idx in ring_atoms[0]
                                        ) or (
                                            end_idx == carbonyl_carbon
                                            and begin_idx in ring_atoms[0]
                                        ):
                                            reactant_has_pyridone = True
                                            break
                                    if reactant_has_pyridone:
                                        break
                                if reactant_has_pyridone:
                                    break

                # Check for pyridone in product
                product_has_pyridone = False
                if checker.check_ring("pyridine", product) and checker.check_fg("Ketone", product):
                    pyridine_indices = checker.get_ring_atom_indices("pyridine", product)
                    ketone_indices = checker.get_fg_atom_indices("Ketone", product)
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        for ring_atoms in pyridine_indices:
                            for ketone_atoms in ketone_indices:
                                carbonyl_carbon = ketone_atoms[0][0]
                                for bond in mol.GetBonds():
                                    begin_idx = bond.GetBeginAtomIdx()
                                    end_idx = bond.GetEndAtomIdx()
                                    if (
                                        begin_idx == carbonyl_carbon and end_idx in ring_atoms[0]
                                    ) or (
                                        end_idx == carbonyl_carbon and begin_idx in ring_atoms[0]
                                    ):
                                        product_has_pyridone = True
                                        break
                                if product_has_pyridone:
                                    break
                            if product_has_pyridone:
                                break

                # In retrosynthesis, if product has pyridone but reactants don't, it means
                # pyridone was created in the forward direction (not preserved)
                if product_has_pyridone and not reactant_has_pyridone:
                    print(f"Pyridone not preserved in reaction (created): {rsmi}")
                    heterocycle_tracking["pyridone"]["preserved"] = False

        # Traverse children
        for child in node.get("children", []):
            track_heterocycles(child, depth + 1)

    # Start tracking
    track_heterocycles(route)

    # Final comprehensive check for heterocycles in the final product
    mol = Chem.MolFromSmiles(final_product_smiles)
    if mol:
        # Final check for triazole
        if not heterocycle_tracking["triazole"]["present_in_final"]:
            triazole_patterns = [
                Chem.MolFromSmarts("c1nn[nH]c1"),  # 1,2,3-triazole
                Chem.MolFromSmarts("c1n[nH]nc1"),  # alternative 1,2,3-triazole
                Chem.MolFromSmarts("c1nnc[nH]1"),  # 1,2,4-triazole
            ]
            for pattern in triazole_patterns:
                if pattern and mol.HasSubstructMatch(pattern):
                    print(f"Found triazole pattern in final product: {final_product_smiles}")
                    heterocycle_tracking["triazole"]["present_in_final"] = True
                    break

        # Final check for pyridone
        if not heterocycle_tracking["pyridone"]["present_in_final"]:
            pyridone_patterns = [
                Chem.MolFromSmarts("c1ccnc(=O)c1"),  # 2-pyridone
                Chem.MolFromSmarts("c1ccn(C)c(=O)c1"),  # N-substituted 2-pyridone
            ]
            for pattern in pyridone_patterns:
                if pattern and mol.HasSubstructMatch(pattern):
                    print(f"Found pyridone pattern in final product: {final_product_smiles}")
                    heterocycle_tracking["pyridone"]["present_in_final"] = True
                    break

    # Return True if either heterocycle is preserved and present in final product
    triazole_result = (
        heterocycle_tracking["triazole"]["preserved"]
        and heterocycle_tracking["triazole"]["present_in_final"]
    )
    pyridone_result = (
        heterocycle_tracking["pyridone"]["preserved"]
        and heterocycle_tracking["pyridone"]["present_in_final"]
    )

    print(
        f"Triazole preserved: {heterocycle_tracking['triazole']['preserved']}, present in final: {heterocycle_tracking['triazole']['present_in_final']}"
    )
    print(
        f"Pyridone preserved: {heterocycle_tracking['pyridone']['preserved']}, present in final: {heterocycle_tracking['pyridone']['present_in_final']}"
    )

    return (
        triazole_result or pyridone_result
    )  # Return True if either heterocycle is preserved and present
