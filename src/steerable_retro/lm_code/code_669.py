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
    This function detects the synthesis of a thiophene-pyrazole biaryl system.
    It checks for a reaction that couples a thiophene-containing compound with a
    pyrazole-containing compound to form a thiophene-pyrazole biaryl system.
    """
    biaryl_synthesis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_synthesis_detected

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains both thiophene and pyrazole rings
                if checker.check_ring("thiophene", product_smiles) and checker.check_ring(
                    "pyrazole", product_smiles
                ):

                    # Check if this is a coupling reaction
                    is_coupling = (
                        checker.check_reaction("Suzuki", rsmi)
                        or checker.check_reaction("Stille", rsmi)
                        or checker.check_reaction("Negishi", rsmi)
                        or checker.check_reaction("Heck", rsmi)
                        or checker.check_reaction("Sonogashira", rsmi)
                        or checker.check_reaction("Kumada", rsmi)
                        or checker.check_reaction("Hiyama-Denmark", rsmi)
                    )

                    if is_coupling:
                        # Check if reactants contain separate rings
                        thiophene_in_reactants = False
                        pyrazole_in_reactants = False

                        for reactant in reactants_smiles:
                            if checker.check_ring("thiophene", reactant):
                                thiophene_in_reactants = True
                            if checker.check_ring("pyrazole", reactant):
                                pyrazole_in_reactants = True

                        # If both rings are present in reactants and it's a coupling reaction,
                        # we've found a thiophene-pyrazole biaryl synthesis
                        if thiophene_in_reactants and pyrazole_in_reactants:
                            print(f"Detected thiophene-pyrazole biaryl synthesis at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            biaryl_synthesis_detected = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Also check if the final product contains the biaryl system
        if node["type"] == "mol" and "smiles" in node and depth == 0:
            mol_smiles = node["smiles"]
            if checker.check_ring("thiophene", mol_smiles) and checker.check_ring(
                "pyrazole", mol_smiles
            ):

                # Check if the rings are connected (biaryl system)
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Get atom indices for both rings
                    thiophene_indices = checker.get_ring_atom_indices("thiophene", mol_smiles)
                    pyrazole_indices = checker.get_ring_atom_indices("pyrazole", mol_smiles)

                    if thiophene_indices and pyrazole_indices:
                        # Check if any atom in thiophene is bonded to any atom in pyrazole
                        for thiophene_atoms in thiophene_indices:
                            for pyrazole_atoms in pyrazole_indices:
                                for t_atom in thiophene_atoms:
                                    for p_atom in pyrazole_atoms:
                                        bond = mol.GetBondBetweenAtoms(t_atom, p_atom)
                                        if bond:
                                            print(
                                                "Detected thiophene-pyrazole biaryl system in final product"
                                            )
                                            biaryl_synthesis_detected = True

        # Continue traversing the route
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return biaryl_synthesis_detected
