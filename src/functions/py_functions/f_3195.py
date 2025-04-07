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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthetic route involves aromatic halogenation.
    """
    aromatic_halogenation = False

    def dfs_traverse(node):
        nonlocal aromatic_halogenation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is an aromatic halogenation reaction using the checker
            if (
                checker.check_reaction("Aromatic fluorination", rsmi)
                or checker.check_reaction("Aromatic chlorination", rsmi)
                or checker.check_reaction("Aromatic bromination", rsmi)
                or checker.check_reaction("Aromatic iodination", rsmi)
            ):

                print(f"Aromatic halogenation reaction detected: {rsmi}")
                aromatic_halogenation = True
                return

            # If not directly identified, check for the transformation
            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has aromatic halide
                if checker.check_fg("Aromatic halide", product_smiles):
                    print(f"Product has aromatic halide: {product_smiles}")

                    # Check if reactants don't have aromatic halide (new formation)
                    reactants_with_aromatic_halide = [
                        r
                        for r in reactants_smiles
                        if checker.check_fg("Aromatic halide", r)
                    ]

                    if len(reactants_with_aromatic_halide) < len(reactants_smiles):
                        # Check for aromatic rings in reactants (something to halogenate)
                        aromatic_rings_in_reactants = False
                        for r_smi in reactants_smiles:
                            if (
                                checker.check_ring("benzene", r_smi)
                                or checker.check_ring("pyridine", r_smi)
                                or checker.check_ring("pyrrole", r_smi)
                                or checker.check_ring("furan", r_smi)
                                or checker.check_ring("thiophene", r_smi)
                                or checker.check_ring("imidazole", r_smi)
                                or checker.check_ring("pyrazole", r_smi)
                                or checker.check_ring("oxazole", r_smi)
                                or checker.check_ring("thiazole", r_smi)
                                or checker.check_ring("indole", r_smi)
                                or checker.check_ring("naphthalene", r_smi)
                            ):
                                aromatic_rings_in_reactants = True
                                break

                        if aromatic_rings_in_reactants:
                            # Check for halogenating agents or halogens in reactants
                            halogen_sources = [
                                "I2",
                                "Br2",
                                "Cl2",
                                "F2",
                                "[I][I]",
                                "[Br][Br]",
                                "[Cl][Cl]",
                                "[F][F]",
                                "NBS",
                                "NCS",
                                "NIS",
                                "ICl",
                                "IBr",
                            ]

                            halogen_present = False
                            for r_smi in reactants_smiles:
                                # Check for common halogenating agents
                                if any(source in r_smi for source in halogen_sources):
                                    halogen_present = True
                                    break

                                # Check for halogen atoms in reactants
                                mol = Chem.MolFromSmiles(r_smi)
                                if mol:
                                    for atom in mol.GetAtoms():
                                        if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                            halogen_present = True
                                            break

                            if halogen_present:
                                print(
                                    f"Aromatic halogenation transformation detected: {rsmi}"
                                )
                                aromatic_halogenation = True
                                return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return aromatic_halogenation
