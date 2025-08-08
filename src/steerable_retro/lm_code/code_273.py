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
    Detects if the synthesis follows a convergent approach where three main fragments
    are combined: a heterocyclic aldehyde, a protected cyclic amine, and a benzyl fragment.
    """
    found_heterocycle = False
    found_protected_amine = False
    found_benzyl_fragment = False

    def dfs_traverse(node):
        nonlocal found_heterocycle, found_protected_amine, found_benzyl_fragment

        if node["type"] == "reaction":
            # Extract reactants and product
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check each reactant for the required fragments
                for reactant_smiles in reactants_smiles:
                    if not reactant_smiles:
                        continue

                    # Check for heterocyclic aldehyde
                    if checker.check_fg("Aldehyde", reactant_smiles) and any(
                        checker.check_ring(ring, reactant_smiles)
                        for ring in [
                            "pyridine",
                            "pyrimidine",
                            "pyrazine",
                            "imidazole",
                            "thiazole",
                            "oxazole",
                            "triazole",
                            "furan",
                            "thiophene",
                            "isoxazole",
                            "isothiazole",
                        ]
                    ):
                        found_heterocycle = True
                        print(f"Found heterocyclic aldehyde: {reactant_smiles}")

                    # Check for protected cyclic amine
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        # Check for various protection groups on cyclic amines
                        has_protection = checker.check_fg(
                            "Boc", reactant_smiles
                        ) or checker.check_fg("Carbamic ester", reactant_smiles)

                        has_cyclic_amine = any(
                            checker.check_ring(ring, reactant_smiles)
                            for ring in [
                                "piperidine",
                                "pyrrolidine",
                                "morpholine",
                                "piperazine",
                                "azetidine",
                                "diazepane",
                                "azepane",
                                "pyrroline",
                            ]
                        )

                        # Check for protection patterns
                        protection_patterns = [
                            "C(=O)OCc1ccccc1",  # Cbz
                            "C(=O)OC(C)(C)C",  # Boc
                            "NC(=O)O[C,c]",  # Carbamate
                            "N(C(=O)O)",  # General carbamate
                            "N(C(=O)OC)",  # N-Carbamate
                        ]

                        has_protection_pattern = False
                        for patt in protection_patterns:
                            patt_mol = Chem.MolFromSmarts(patt)
                            if patt_mol and reactant_mol.HasSubstructMatch(patt_mol):
                                has_protection_pattern = True
                                break

                        # Check for nitrogen in a ring
                        nitrogen_in_ring = False
                        for atom in reactant_mol.GetAtoms():
                            if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                nitrogen_in_ring = True
                                break

                        if (has_protection or has_protection_pattern) and (
                            has_cyclic_amine or nitrogen_in_ring
                        ):
                            found_protected_amine = True
                            print(f"Found protected cyclic amine: {reactant_smiles}")

                    # Check for benzyl fragment
                    if reactant_mol:
                        # Look for benzyl group with potential substituents
                        benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[!#1]")
                        if reactant_mol.HasSubstructMatch(benzyl_pattern):
                            # Verify it's a relevant benzyl fragment, not just any benzyl group
                            # Check if it has a functional group attached or is a halide
                            if (
                                checker.check_fg("Primary halide", reactant_smiles)
                                or checker.check_fg("Secondary halide", reactant_smiles)
                                or "Br" in reactant_smiles
                                or "Cl" in reactant_smiles
                                or "I" in reactant_smiles
                            ):
                                found_benzyl_fragment = True
                                print(f"Found benzyl fragment: {reactant_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if all three fragments are found
    result = found_heterocycle and found_protected_amine and found_benzyl_fragment
    print(f"Found heterocyclic aldehyde: {found_heterocycle}")
    print(f"Found protected cyclic amine: {found_protected_amine}")
    print(f"Found benzyl fragment: {found_benzyl_fragment}")
    print(f"Final result: {result}")

    return result
