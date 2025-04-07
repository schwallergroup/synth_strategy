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
    Detects a strategy involving the coupling of a nitrile-containing fragment
    in the late stages of the synthesis.
    """
    found_nitrile_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_coupling

        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Focus on late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if any reactant contains a nitrile group
                reactant_smiles_list = [r for r in reactants_str.split(".") if r]
                product_smiles = product_str if product_str else ""

                # Find nitrile-containing reactants
                nitrile_reactants = []
                for r_smiles in reactant_smiles_list:
                    if checker.check_fg("Nitrile", r_smiles):
                        nitrile_reactants.append(r_smiles)
                        print(f"Found nitrile-containing reactant: {r_smiles}")

                # Check if product also contains nitrile (preservation)
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)
                print(f"Product has nitrile: {product_has_nitrile}")

                if nitrile_reactants and product_has_nitrile:
                    # Check if this is a coupling reaction
                    is_coupling = False

                    # Check for common coupling reaction types
                    coupling_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic esters",
                        "Suzuki coupling with boronic acids OTf",
                        "Suzuki coupling with boronic esters OTf",
                        "Suzuki coupling with sulfonic esters",
                        "Negishi coupling",
                        "Stille reaction_aryl",
                        "Stille reaction_vinyl",
                        "Stille reaction_benzyl",
                        "Stille reaction_allyl",
                        "Stille reaction_aryl OTf",
                        "Stille reaction_vinyl OTf",
                        "Stille reaction_benzyl OTf",
                        "Stille reaction_allyl OTf",
                        "Sonogashira alkyne_aryl halide",
                        "Sonogashira acetylene_aryl halide",
                        "Sonogashira alkyne_aryl OTf",
                        "Sonogashira acetylene_aryl OTf",
                        "Heck terminal vinyl",
                        "Heck_non-terminal_vinyl",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "Ullmann condensation",
                        "Hiyama-Denmark Coupling",
                        "Kumada cross-coupling",
                        "Aryllithium cross-coupling",
                        "decarboxylative_coupling",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "N-arylation_heterocycles",
                        "Ullmann-Goldberg Substitution amine",
                        "Ullmann-Goldberg Substitution thiol",
                        "Ullmann-Goldberg Substitution aryl alcohol",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_coupling = True
                            print(f"Identified coupling reaction: {rxn_type}")
                            break

                    # Check for significant fragments (alternative approach)
                    if not is_coupling:
                        print("Checking for significant fragments...")
                        significant_fragments = 0
                        for r_smiles in reactant_smiles_list:
                            r_mol = Chem.MolFromSmiles(r_smiles)
                            if (
                                r_mol and r_mol.GetNumHeavyAtoms() >= 5
                            ):  # Consider fragments with at least 5 heavy atoms
                                significant_fragments += 1
                                print(
                                    f"Found significant fragment: {r_smiles} with {r_mol.GetNumHeavyAtoms()} heavy atoms"
                                )

                        if significant_fragments >= 2 or (
                            len(reactant_smiles_list) >= 2
                            and significant_fragments >= 1
                        ):
                            is_coupling = True
                            print(f"Identified coupling based on significant fragments")

                    # Check for C-C bond formation (another alternative approach)
                    if not is_coupling and len(reactant_smiles_list) >= 2:
                        print("Checking for C-C bond formation...")
                        try:
                            # Check if there are aromatic rings in both reactants and product
                            aromatic_reactants = 0
                            for r_smiles in reactant_smiles_list:
                                r_mol = Chem.MolFromSmiles(r_smiles)
                                if r_mol:
                                    ring_info = r_mol.GetRingInfo()
                                    if ring_info.NumRings() > 0:
                                        aromatic_reactants += 1

                            product_mol = Chem.MolFromSmiles(product_smiles)
                            if product_mol and aromatic_reactants >= 2:
                                product_rings = product_mol.GetRingInfo().NumRings()
                                # If product has rings and multiple reactants have rings, likely a coupling
                                if product_rings > 0:
                                    is_coupling = True
                                    print(f"Identified coupling based on ring analysis")
                        except Exception as e:
                            print(f"Error in C-C bond analysis: {e}")

                    # Check for common rings in the product that might indicate coupling
                    if not is_coupling:
                        print("Checking for biaryl or heterocycle formation...")
                        common_coupling_rings = [
                            "benzene",
                            "pyridine",
                            "pyrimidine",
                            "pyrazine",
                            "pyridazine",
                            "thiophene",
                            "furan",
                            "pyrrole",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                        ]

                        for ring_name in common_coupling_rings:
                            if checker.check_ring(ring_name, product_smiles):
                                # Check if this ring is newly formed
                                ring_in_reactants = False
                                for r_smiles in reactant_smiles_list:
                                    if checker.check_ring(ring_name, r_smiles):
                                        ring_in_reactants = True
                                        break

                                if not ring_in_reactants:
                                    is_coupling = True
                                    print(
                                        f"Identified coupling based on formation of {ring_name} ring"
                                    )
                                    break

                    if is_coupling:
                        found_nitrile_coupling = True
                        print(
                            f"Found nitrile-containing fragment coupling at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_nitrile_coupling}")

    return found_nitrile_coupling
