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
    This function detects a synthesis strategy involving the coupling of a pyrazole-containing
    fragment with a difluoromethylbenzene-containing fragment.
    """
    has_pyrazole_reactant = False
    has_difluoro_reactant = False
    has_coupling = False

    # More general patterns for difluorinated benzenes
    difluoromethyl_pattern = Chem.MolFromSmarts("c-C(F)(F)")  # Difluoromethyl directly on benzene
    difluorobenzene_pattern = Chem.MolFromSmarts("c1c(F)cc(F)cc1")  # 1,3-difluorobenzene
    difluorobenzene_pattern2 = Chem.MolFromSmarts("c1cc(F)c(F)cc1")  # 1,2-difluorobenzene
    difluorobenzene_pattern3 = Chem.MolFromSmarts("c1cc(F)cc(F)c1")  # 1,4-difluorobenzene

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole_reactant, has_difluoro_reactant, has_coupling

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for pyrazole
            if checker.check_ring("pyrazole", mol_smiles):
                print(f"Found pyrazole in molecule: {mol_smiles}")
                has_pyrazole_reactant = True

            # Check for difluorinated benzene structures
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                if (
                    mol.HasSubstructMatch(difluoromethyl_pattern)
                    or mol.HasSubstructMatch(difluorobenzene_pattern)
                    or mol.HasSubstructMatch(difluorobenzene_pattern2)
                    or mol.HasSubstructMatch(difluorobenzene_pattern3)
                ):
                    print(f"Found difluorinated benzene in molecule: {mol_smiles}")
                    has_difluoro_reactant = True

                # Additional check for any benzene with two fluorines
                fluorine_count = 0
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F":
                        fluorine_count += 1

                if (
                    fluorine_count >= 2 and mol_smiles.find("c") >= 0
                ):  # Has at least 2 fluorines and aromatic carbon
                    print(
                        f"Found molecule with multiple fluorines and aromatic system: {mol_smiles}"
                    )
                    has_difluoro_reactant = True

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain the fragments separately
                has_pyrazole_in_reactants = False
                has_difluoro_in_reactants = False

                for reactant_smiles in reactants_smiles:
                    if checker.check_ring("pyrazole", reactant_smiles):
                        has_pyrazole_in_reactants = True

                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        if (
                            reactant_mol.HasSubstructMatch(difluoromethyl_pattern)
                            or reactant_mol.HasSubstructMatch(difluorobenzene_pattern)
                            or reactant_mol.HasSubstructMatch(difluorobenzene_pattern2)
                            or reactant_mol.HasSubstructMatch(difluorobenzene_pattern3)
                        ):
                            has_difluoro_in_reactants = True

                        # Count fluorines
                        fluorine_count = 0
                        for atom in reactant_mol.GetAtoms():
                            if atom.GetSymbol() == "F":
                                fluorine_count += 1

                        if fluorine_count >= 2 and reactant_smiles.find("c") >= 0:
                            has_difluoro_in_reactants = True

                # Check if product contains both fragments
                product_has_pyrazole = checker.check_ring("pyrazole", product_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)
                product_has_difluoro = False

                if product_mol:
                    if (
                        product_mol.HasSubstructMatch(difluoromethyl_pattern)
                        or product_mol.HasSubstructMatch(difluorobenzene_pattern)
                        or product_mol.HasSubstructMatch(difluorobenzene_pattern2)
                        or product_mol.HasSubstructMatch(difluorobenzene_pattern3)
                    ):
                        product_has_difluoro = True

                    # Count fluorines
                    fluorine_count = 0
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() == "F":
                            fluorine_count += 1

                    if fluorine_count >= 2 and product_smiles.find("c") >= 0:
                        product_has_difluoro = True

                # Check for coupling reactions
                is_coupling_reaction = (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                )

                # Verify coupling: either explicit coupling reaction or fragments joined
                if is_coupling_reaction:
                    print(f"Found explicit coupling reaction at depth {depth}: {rsmi}")
                    has_coupling = True
                elif (
                    has_pyrazole_in_reactants
                    and has_difluoro_in_reactants
                    and product_has_pyrazole
                    and product_has_difluoro
                ):
                    print(
                        f"Found reaction joining pyrazole and difluoro fragments at depth {depth}: {rsmi}"
                    )
                    has_coupling = True
                elif len(reactants_smiles) >= 2 and product_has_pyrazole and product_has_difluoro:
                    # Check if one reactant has pyrazole and another has difluoro
                    pyrazole_reactant = None
                    difluoro_reactant = None

                    for reactant_smiles in reactants_smiles:
                        if checker.check_ring("pyrazole", reactant_smiles):
                            pyrazole_reactant = reactant_smiles

                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol:
                            has_difluoro = False
                            if (
                                reactant_mol.HasSubstructMatch(difluoromethyl_pattern)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern2)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern3)
                            ):
                                has_difluoro = True

                            # Count fluorines
                            fluorine_count = 0
                            for atom in reactant_mol.GetAtoms():
                                if atom.GetSymbol() == "F":
                                    fluorine_count += 1

                            if fluorine_count >= 2 and reactant_smiles.find("c") >= 0:
                                has_difluoro = True

                            if has_difluoro:
                                difluoro_reactant = reactant_smiles

                    if (
                        pyrazole_reactant
                        and difluoro_reactant
                        and pyrazole_reactant != difluoro_reactant
                    ):
                        print(
                            f"Found reaction combining separate pyrazole and difluoro fragments at depth {depth}: {rsmi}"
                        )
                        has_coupling = True

                # Also check for amide coupling specifically
                if "C(=O)" in rsmi and "N" in rsmi and len(reactants_smiles) >= 2:
                    # One reactant has pyrazole, another has difluoro, product has both
                    has_pyrazole_only = False
                    has_difluoro_only = False

                    for reactant_smiles in reactants_smiles:
                        has_pyrazole = checker.check_ring("pyrazole", reactant_smiles)

                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        has_difluoro = False
                        if reactant_mol:
                            if (
                                reactant_mol.HasSubstructMatch(difluoromethyl_pattern)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern2)
                                or reactant_mol.HasSubstructMatch(difluorobenzene_pattern3)
                            ):
                                has_difluoro = True

                            # Count fluorines
                            fluorine_count = 0
                            for atom in reactant_mol.GetAtoms():
                                if atom.GetSymbol() == "F":
                                    fluorine_count += 1

                            if fluorine_count >= 2 and reactant_smiles.find("c") >= 0:
                                has_difluoro = True

                        if has_pyrazole and not has_difluoro:
                            has_pyrazole_only = True
                        elif has_difluoro and not has_pyrazole:
                            has_difluoro_only = True

                    if (
                        has_pyrazole_only
                        and has_difluoro_only
                        and product_has_pyrazole
                        and product_has_difluoro
                    ):
                        print(
                            f"Found amide coupling joining pyrazole and difluoro fragments at depth {depth}: {rsmi}"
                        )
                        has_coupling = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Results: pyrazole={has_pyrazole_reactant}, difluoro={has_difluoro_reactant}, coupling={has_coupling}"
    )
    return has_pyrazole_reactant and has_difluoro_reactant and has_coupling
