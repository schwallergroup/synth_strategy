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
    This function detects if the synthesis follows a convergent approach with
    multiple fragments being joined in the second half of the synthesis.
    """
    is_convergent = False
    max_depth = 0

    # First pass to determine the maximum depth of the synthesis
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")

    # Second pass to detect convergent synthesis
    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a reaction with multiple reactants
            if len(reactants) > 1 and all(r.strip() for r in reactants):
                # Check if this is in the second half of synthesis
                if depth <= (max_depth // 2 + 1):
                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    # Check if this is a coupling reaction
                    coupling_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic esters",
                        "Suzuki coupling with boronic acids OTf",
                        "Suzuki coupling with boronic esters OTf",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "Buchwald-Hartwig",
                        "Sonogashira acetylene_aryl halide",
                        "Sonogashira alkyne_aryl halide",
                        "Negishi coupling",
                        "Stille reaction_aryl",
                        "Stille reaction_vinyl",
                        "Heck terminal vinyl",
                        "Ullmann-Goldberg Substitution amine",
                        "Goldberg coupling",
                    ]

                    is_coupling = False
                    for rxn in coupling_reactions:
                        if checker.check_reaction(rxn, rsmi):
                            print(f"Detected coupling reaction: {rxn}")
                            is_coupling = True
                            break

                    # If not identified as a coupling reaction, check for other reactions
                    if not is_coupling:
                        other_reactions = [
                            "Wittig reaction with triphenylphosphorane",
                            "Wittig with Phosphonium",
                            "Grignard from aldehyde to alcohol",
                            "Grignard from ketone to alcohol",
                            "Mitsunobu aryl ether",
                            "Diels-Alder",
                            "Michael addition",
                            "Aldol condensation",
                            "Ugi reaction",
                        ]

                        for rxn in other_reactions:
                            if checker.check_reaction(rxn, rsmi):
                                print(f"Detected other convergent reaction: {rxn}")
                                is_coupling = True
                                break

                    # If still not identified, check for C-C bond formation between aromatic rings
                    if not is_coupling:
                        product = rsmi.split(">")[-1]
                        product_mol = Chem.MolFromSmiles(product)

                        # Check if the product contains a biaryl system
                        if (
                            product_mol
                            and checker.check_fg("Aromatic halide", reactants[0])
                            and any("B" in r or "Sn" in r or "Zn" in r for r in reactants)
                        ):
                            print(
                                "Detected potential coupling reaction based on reactants and product structure"
                            )
                            is_coupling = True

                    print(f"Is coupling reaction: {is_coupling}")

                    if is_coupling:
                        # Count complex reactants (significant fragments)
                        complex_reactants = 0
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                num_heavy_atoms = mol.GetNumHeavyAtoms()
                                num_rings = rdMolDescriptors.CalcNumRings(mol)
                                print(
                                    f"Reactant: {reactant}, Heavy atoms: {num_heavy_atoms}, Rings: {num_rings}"
                                )

                                # Consider reactants with significant structure as complex
                                if num_heavy_atoms >= 6 or num_rings > 0:
                                    complex_reactants += 1

                        print(f"Complex reactants: {complex_reactants}")
                        if complex_reactants >= 2:
                            print(f"Convergent synthesis detected at depth {depth}")
                            is_convergent = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return is_convergent
