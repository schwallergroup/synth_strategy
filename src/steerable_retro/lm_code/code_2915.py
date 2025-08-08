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
    Detects the presence of a four-carbon linker between indole and piperazine-pyrimidine fragments.
    """
    found_four_carbon_linker = False

    def dfs_traverse(node, depth=0):
        nonlocal found_four_carbon_linker

        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]

            # Check if molecule contains all required rings
            has_indole = checker.check_ring("indole", mol_smiles)
            has_piperazine = checker.check_ring("piperazine", mol_smiles)
            has_pyrimidine = checker.check_ring("pyrimidine", mol_smiles)

            # If all required rings are present, check for the four-carbon linker
            if has_indole and has_piperazine and has_pyrimidine:
                print(f"Found molecule with indole, piperazine, and pyrimidine: {mol_smiles}")

                # Create RDKit molecule
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    try:
                        # Get atom indices for each ring
                        indole_indices_list = checker.get_ring_atom_indices("indole", mol_smiles)
                        piperazine_indices_list = checker.get_ring_atom_indices(
                            "piperazine", mol_smiles
                        )

                        # Create sets of atoms for each ring
                        indole_atoms = set()
                        if indole_indices_list:
                            for indices_tuple in indole_indices_list:
                                if isinstance(indices_tuple[0], tuple):
                                    for idx in indices_tuple[0]:
                                        indole_atoms.add(idx)
                                else:
                                    for idx in indices_tuple:
                                        indole_atoms.add(idx)

                        piperazine_atoms = set()
                        if piperazine_indices_list:
                            for indices_tuple in piperazine_indices_list:
                                if isinstance(indices_tuple[0], tuple):
                                    for idx in indices_tuple[0]:
                                        piperazine_atoms.add(idx)
                                else:
                                    for idx in indices_tuple:
                                        piperazine_atoms.add(idx)

                        # If we have atoms for both rings, check for paths
                        if indole_atoms and piperazine_atoms:
                            # Check for paths between indole and piperazine
                            for indole_atom in indole_atoms:
                                for piperazine_atom in piperazine_atoms:
                                    # Find all paths of length 5 (4 bonds) between the atoms
                                    paths = Chem.FindAllPathsOfLengthN(
                                        mol, 5, indole_atom, piperazine_atom
                                    )

                                    for path in paths:
                                        # Check if the path is a 4-carbon chain (5 atoms including start and end)
                                        if len(path) == 5:
                                            # Verify the middle atoms are carbons (atomic number 6)
                                            middle_atoms = path[1:-1]
                                            if all(
                                                mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
                                                for idx in middle_atoms
                                            ):
                                                print(
                                                    f"Found four-carbon linker between indole and piperazine in: {mol_smiles}"
                                                )
                                                print(f"Path: {path}")
                                                found_four_carbon_linker = True
                                                return  # Exit early once found
                    except Exception as e:
                        print(f"Error analyzing molecule structure: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if the product has all required rings
                product_has_indole = checker.check_ring("indole", product)
                product_has_piperazine = checker.check_ring("piperazine", product)
                product_has_pyrimidine = checker.check_ring("pyrimidine", product)

                if product_has_indole and product_has_piperazine and product_has_pyrimidine:
                    # Create RDKit molecule for product
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol:
                        try:
                            # Check for 4-carbon linker in product
                            product_has_linker = False

                            # Get atom indices for each ring in product
                            indole_indices_list = checker.get_ring_atom_indices("indole", product)
                            piperazine_indices_list = checker.get_ring_atom_indices(
                                "piperazine", product
                            )

                            # Create sets of atoms for each ring
                            indole_atoms = set()
                            if indole_indices_list:
                                for indices_tuple in indole_indices_list:
                                    if isinstance(indices_tuple[0], tuple):
                                        for idx in indices_tuple[0]:
                                            indole_atoms.add(idx)
                                    else:
                                        for idx in indices_tuple:
                                            indole_atoms.add(idx)

                            piperazine_atoms = set()
                            if piperazine_indices_list:
                                for indices_tuple in piperazine_indices_list:
                                    if isinstance(indices_tuple[0], tuple):
                                        for idx in indices_tuple[0]:
                                            piperazine_atoms.add(idx)
                                    else:
                                        for idx in indices_tuple:
                                            piperazine_atoms.add(idx)

                            # If we have atoms for both rings, check for paths
                            if indole_atoms and piperazine_atoms:
                                # Check for paths between indole and piperazine in product
                                for indole_atom in indole_atoms:
                                    for piperazine_atom in piperazine_atoms:
                                        paths = Chem.FindAllPathsOfLengthN(
                                            product_mol, 5, indole_atom, piperazine_atom
                                        )

                                        for path in paths:
                                            if len(path) == 5:
                                                middle_atoms = path[1:-1]
                                                if all(
                                                    product_mol.GetAtomWithIdx(idx).GetAtomicNum()
                                                    == 6
                                                    for idx in middle_atoms
                                                ):
                                                    product_has_linker = True
                                                    break

                                    if product_has_linker:
                                        break

                                if product_has_linker:
                                    # Check if any reactant already has the complete structure
                                    reactants_with_linker = 0
                                    for reactant in reactants:
                                        if (
                                            checker.check_ring("indole", reactant)
                                            and checker.check_ring("piperazine", reactant)
                                            and checker.check_ring("pyrimidine", reactant)
                                        ):

                                            reactant_mol = Chem.MolFromSmiles(reactant)
                                            if reactant_mol:
                                                try:
                                                    # Get atom indices for each ring in reactant
                                                    r_indole_indices = (
                                                        checker.get_ring_atom_indices(
                                                            "indole", reactant
                                                        )
                                                    )
                                                    r_piperazine_indices = (
                                                        checker.get_ring_atom_indices(
                                                            "piperazine", reactant
                                                        )
                                                    )

                                                    # Create sets of atoms for each ring
                                                    r_indole_atoms = set()
                                                    if r_indole_indices:
                                                        for indices_tuple in r_indole_indices:
                                                            if isinstance(indices_tuple[0], tuple):
                                                                for idx in indices_tuple[0]:
                                                                    r_indole_atoms.add(idx)
                                                            else:
                                                                for idx in indices_tuple:
                                                                    r_indole_atoms.add(idx)

                                                    r_piperazine_atoms = set()
                                                    if r_piperazine_indices:
                                                        for indices_tuple in r_piperazine_indices:
                                                            if isinstance(indices_tuple[0], tuple):
                                                                for idx in indices_tuple[0]:
                                                                    r_piperazine_atoms.add(idx)
                                                            else:
                                                                for idx in indices_tuple:
                                                                    r_piperazine_atoms.add(idx)

                                                    # If we have atoms for both rings, check for paths
                                                    if r_indole_atoms and r_piperazine_atoms:
                                                        for r_indole_atom in r_indole_atoms:
                                                            for (
                                                                r_piperazine_atom
                                                            ) in r_piperazine_atoms:
                                                                paths = Chem.FindAllPathsOfLengthN(
                                                                    reactant_mol,
                                                                    5,
                                                                    r_indole_atom,
                                                                    r_piperazine_atom,
                                                                )

                                                                for path in paths:
                                                                    if len(path) == 5:
                                                                        middle_atoms = path[1:-1]
                                                                        if all(
                                                                            reactant_mol.GetAtomWithIdx(
                                                                                idx
                                                                            ).GetAtomicNum()
                                                                            == 6
                                                                            for idx in middle_atoms
                                                                        ):
                                                                            reactants_with_linker += (
                                                                                1
                                                                            )
                                                                            break
                                                except Exception as e:
                                                    print(
                                                        f"Error analyzing reactant structure: {e}"
                                                    )

                                    if reactants_with_linker == 0:
                                        print(
                                            f"Found reaction forming four-carbon linker between indole and piperazine-pyrimidine: {rsmi}"
                                        )
                                        found_four_carbon_linker = True
                                        return  # Exit early once found
                        except Exception as e:
                            print(f"Error analyzing product structure: {e}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)
            if found_four_carbon_linker:
                return  # Exit early if pattern found in any child

    # Start traversal
    dfs_traverse(route)

    return found_four_carbon_linker
