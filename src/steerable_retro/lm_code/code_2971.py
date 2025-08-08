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
    This function detects a specific functional group interconversion sequence:
    nitro → amine → isothiocyanate → heterocycle
    """
    # Initialize flags for each transformation
    nitro_reduction = False
    isothiocyanate_formation = False
    heterocycle_formation = False

    # Track atom mappings to ensure we're following the same functional group
    nitro_atom_maps = set()
    amine_atom_maps = set()
    isothiocyanate_atom_maps = set()

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction, isothiocyanate_formation, heterocycle_formation
        nonlocal nitro_atom_maps, amine_atom_maps, isothiocyanate_atom_maps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for nitro reduction (depth 5)
                if depth == 5:
                    # Check if this is a nitro reduction reaction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print("Detected nitro reduction reaction at depth 5")

                        # Verify nitro group in reactant and amine in product
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant):
                                print(f"Found nitro group in reactant: {reactant}")

                                # Get the product and check for amine
                                if checker.check_fg("Primary amine", product):
                                    print(f"Found primary amine in product: {product}")
                                    nitro_reduction = True

                                    # Extract atom mappings for tracking
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    product_mol = Chem.MolFromSmiles(product)

                                    if reactant_mol and product_mol:
                                        for atom in reactant_mol.GetAtoms():
                                            if atom.GetSymbol() == "N" and atom.HasProp(
                                                "molAtomMapNumber"
                                            ):
                                                nitro_atom_maps.add(
                                                    atom.GetProp("molAtomMapNumber")
                                                )

                                        for atom in product_mol.GetAtoms():
                                            if atom.GetSymbol() == "N" and atom.HasProp(
                                                "molAtomMapNumber"
                                            ):
                                                amine_atom_maps.add(
                                                    atom.GetProp("molAtomMapNumber")
                                                )

                                    print(f"Nitro atom maps: {nitro_atom_maps}")
                                    print(f"Amine atom maps: {amine_atom_maps}")
                    else:
                        # Fallback: check for nitro to amine conversion without specific reaction
                        nitro_in_reactants = any(
                            checker.check_fg("Nitro group", r) for r in reactants
                        )
                        amine_in_product = checker.check_fg("Primary amine", product)

                        if nitro_in_reactants and amine_in_product:
                            print("Detected nitro to amine conversion at depth 5 (fallback)")
                            nitro_reduction = True

                # Check for isothiocyanate formation (depth 3)
                if depth == 3:
                    # Check if this is an isothiocyanate formation reaction
                    if checker.check_reaction("Amine and thiophosgene to isothiocyanate", rsmi):
                        print("Detected isothiocyanate formation reaction at depth 3")

                        # Verify amine in reactants and isothiocyanate in product
                        amine_found = False
                        for reactant in reactants:
                            if checker.check_fg("Primary amine", reactant):
                                print(f"Found primary amine in reactant: {reactant}")
                                amine_found = True

                        if amine_found and checker.check_fg("Isothiocyanate", product):
                            print(f"Found isothiocyanate in product: {product}")
                            isothiocyanate_formation = True
                    else:
                        # Fallback: check for amine to isothiocyanate conversion without specific reaction
                        amine_in_reactants = any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        )
                        isothiocyanate_in_product = checker.check_fg("Isothiocyanate", product)

                        if amine_in_reactants and isothiocyanate_in_product:
                            print(
                                "Detected amine to isothiocyanate conversion at depth 3 (fallback)"
                            )
                            isothiocyanate_formation = True

                # Check for heterocycle formation from isothiocyanate (depth 1)
                if depth == 1:
                    # Check for heterocycle formation reactions
                    heterocycle_reactions = [
                        "Benzimidazole formation from aldehyde",
                        "Benzimidazole formation from acyl halide",
                        "Benzimidazole formation from ester/carboxylic acid",
                        "Benzothiazole formation from aldehyde",
                        "Benzothiazole formation from acyl halide",
                        "Benzothiazole formation from ester/carboxylic acid",
                        "Benzoxazole formation from aldehyde",
                        "Benzoxazole formation from acyl halide",
                        "Benzoxazole formation from ester/carboxylic acid",
                        "Benzoxazole formation (intramolecular)",
                    ]

                    reaction_detected = False
                    for rxn_type in heterocycle_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected heterocycle formation reaction: {rxn_type}")
                            reaction_detected = True
                            break

                    # Verify isothiocyanate in reactants
                    isothiocyanate_found = any(
                        checker.check_fg("Isothiocyanate", r) for r in reactants
                    )
                    if isothiocyanate_found:
                        print("Found isothiocyanate in reactants")

                    # Check for heterocycle in product
                    heterocycle_rings = [
                        "benzimidazole",
                        "benzothiazole",
                        "benzoxazole",
                        "thiazole",
                        "oxazole",
                        "imidazole",
                    ]

                    heterocycle_found = False
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product):
                            print(f"Found {ring} in product: {product}")
                            heterocycle_found = True
                            break

                    if (reaction_detected or heterocycle_found) and isothiocyanate_found:
                        heterocycle_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Nitro reduction: {nitro_reduction}")
    print(f"Isothiocyanate formation: {isothiocyanate_formation}")
    print(f"Heterocycle formation: {heterocycle_formation}")

    # Return True only if all transformations are detected in the correct sequence
    return nitro_reduction and isothiocyanate_formation and heterocycle_formation
