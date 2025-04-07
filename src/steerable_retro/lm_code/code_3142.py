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
    This function detects if the synthesis route involves halogenation
    (fluorination, chlorination, bromination, or iodination) in the final step.
    """
    final_step_halogenation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_halogenation

        # Print node information for debugging
        if node["type"] == "mol":
            print(f"Depth {depth}, Molecule: {node['smiles']}")
        elif node["type"] == "reaction":
            print(f"Depth {depth}, Reaction: {node['metadata'].get('rsmi', 'No RSMI')}")

        # Check if this is a reaction at depth 1 (final synthetic step)
        if node["type"] == "reaction" and depth == 1:
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing final step reaction: {rsmi}")

                # Check if this is a halogenation reaction using the checker functions
                halogenation_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Fluorination",
                    "Iodination",
                    "Bromination",
                ]

                for reaction_type in halogenation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} in final step")
                        final_step_halogenation = True
                        return

                # If no specific halogenation reaction was detected, check for halogen addition
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Create RDKit molecules
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                if not product_mol or any(r is None for r in reactant_mols):
                    print("Failed to parse some molecules")
                    return

                # Check for common brominating agents like NBS
                if any("O=C1CCC(=O)N1Br" in r for r in reactants_smiles):
                    print("Detected N-bromosuccinimide (NBS) as a reactant, likely a bromination")
                    final_step_halogenation = True
                    return

                # Check for other common halogenation reagents
                halogenation_reagents = {
                    "F": ["F2", "XeF2", "CF3OF", "CH3COOF"],
                    "Cl": ["Cl2", "NaOCl", "Ca(OCl)2", "SO2Cl2", "PCl5"],
                    "Br": ["Br2", "NBS", "CBr4", "PBr3"],
                    "I": ["I2", "NIS", "ICl"],
                }

                for halogen, reagents in halogenation_reagents.items():
                    for reagent in reagents:
                        for reactant in reactants_smiles:
                            if reagent in reactant:
                                print(f"Detected {halogen}ation reagent: {reagent}")
                                final_step_halogenation = True
                                return

                # Count halogen atoms in product and reactants
                halogen_symbols = ["F", "Cl", "Br", "I"]

                # Count halogens in product
                product_halogens = {}
                for atom in product_mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    if symbol in halogen_symbols:
                        product_halogens[symbol] = product_halogens.get(symbol, 0) + 1

                # Count halogens in reactants
                reactants_halogens = {}
                for reactant in reactant_mols:
                    for atom in reactant.GetAtoms():
                        symbol = atom.GetSymbol()
                        if symbol in halogen_symbols:
                            reactants_halogens[symbol] = reactants_halogens.get(symbol, 0) + 1

                print(f"Product halogens: {product_halogens}")
                print(f"Reactants halogens: {reactants_halogens}")

                # Check if any halogen was added or increased in the product
                for halogen in halogen_symbols:
                    product_count = product_halogens.get(halogen, 0)
                    reactant_count = reactants_halogens.get(halogen, 0)

                    if product_count > reactant_count:
                        print(
                            f"Detected {halogen} addition in final step: {reactant_count} → {product_count}"
                        )
                        final_step_halogenation = True
                        return

                # Check for halogen transfer reactions (like NBS bromination)
                # Look at atom mapping to see if bromine moved to a new position
                if "Br" in rsmi:
                    # Check if the reaction involves a bromine atom with atom mapping
                    # This is a simplification - in a real implementation, we would trace the atom mapping
                    if "[Br:" in rsmi and "Br]" in product_smiles:
                        print("Detected bromine transfer in final step based on atom mapping")
                        final_step_halogenation = True
                        return

                # Also check for halogen functional groups as a backup
                halogen_fgs = [
                    "Aromatic halide",
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Alkenyl halide",
                    "Haloalkyne",
                ]

                product_halogen_fgs = set()
                for fg in halogen_fgs:
                    if checker.check_fg(fg, product_smiles):
                        product_halogen_fgs.add(fg)

                reactants_halogen_fgs = set()
                for reactant in reactants_smiles:
                    for fg in halogen_fgs:
                        if checker.check_fg(fg, reactant):
                            reactants_halogen_fgs.add(fg)

                print(f"Product halogen FGs: {product_halogen_fgs}")
                print(f"Reactants halogen FGs: {reactants_halogen_fgs}")

                # If product has halogen FGs not in reactants, it's likely halogenation
                if product_halogen_fgs - reactants_halogen_fgs:
                    print(
                        f"Detected new halogen functional groups: {product_halogen_fgs - reactants_halogen_fgs}"
                    )
                    final_step_halogenation = True
                    return

                # Special case for NBS bromination where the bromine is transferred
                # Look for patterns in the reaction SMILES that indicate bromination
                if "O=C1CCC(=O)N1[Br:" in rsmi and "[Br:" in product_smiles:
                    print("Detected NBS bromination pattern with atom mapping")
                    final_step_halogenation = True
                    return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return final_step_halogenation
