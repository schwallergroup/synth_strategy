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
    This function detects late-stage N-alkylation,
    specifically N-methylation of a heterocyclic nitrogen.
    """
    n_alkylation_detected = False

    # List of common heterocyclic rings to check
    heterocyclic_rings = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "benzimidazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "pyrrolidine",
    ]

    # List of N-alkylation reaction types
    methylation_reactions = [
        "N-methylation",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "DMS Amine methylation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction" and depth <= 3:  # Late stage (low depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is an N-alkylation reaction
                is_methylation = False
                for rxn_type in methylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected methylation reaction: {rxn_type}")
                        is_methylation = True
                        break

                # Check for general N-alkylation if specific methylation not found
                if not is_methylation:
                    is_alkylation = checker.check_reaction("N-alkylation", rsmi)
                    if is_alkylation:
                        print(f"Detected general N-alkylation reaction")
                        is_methylation = True

                # If we still haven't found a match, check for Buchwald-Hartwig/N-arylation
                if not is_methylation:
                    if (
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                    ):
                        print(f"Detected N-arylation reaction")
                        is_methylation = True

                if is_methylation or "N-alkylation" in rsmi:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if the product contains a heterocyclic ring
                    product_heterocycles = []
                    for ring in heterocyclic_rings:
                        if checker.check_ring(ring, product):
                            product_heterocycles.append(ring)

                    if product_heterocycles:
                        print(f"Product contains heterocycles: {', '.join(product_heterocycles)}")

                        # Check if any reactant has a heterocyclic ring with N-H
                        for reactant in reactants:
                            reactant_heterocycles = []
                            for ring in heterocyclic_rings:
                                if checker.check_ring(ring, reactant):
                                    reactant_heterocycles.append(ring)

                            if reactant_heterocycles:
                                print(
                                    f"Reactant contains heterocycles: {', '.join(reactant_heterocycles)}"
                                )

                                # Check for N-H in heterocyclic rings
                                try:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        # Check for heterocyclic N-H
                                        has_heterocyclic_nh = False
                                        for atom in reactant_mol.GetAtoms():
                                            if (
                                                atom.GetSymbol() == "N"
                                                and atom.IsInRing()
                                                and atom.GetTotalNumHs() > 0
                                            ):
                                                print(f"Found heterocyclic N-H in reactant")
                                                has_heterocyclic_nh = True
                                                break

                                        # Also check using functional group detection
                                        if (
                                            has_heterocyclic_nh
                                            or "[nH]" in reactant
                                            or checker.check_fg("Secondary amine", reactant)
                                            or checker.check_fg("Primary amine", reactant)
                                            or checker.check_fg("Tertiary amine", reactant)
                                        ):
                                            print(
                                                f"Late-stage N-alkylation detected at depth {depth}"
                                            )
                                            print(f"Reaction SMILES: {rsmi}")
                                            n_alkylation_detected = True
                                            return
                                except Exception as e:
                                    print(f"Error processing reactant: {e}")

                # Direct check for N-alkylation by examining the reaction
                if not n_alkylation_detected and depth <= 3:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if product has N-CH3 that wasn't in reactants
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            for atom in product_mol.GetAtoms():
                                if atom.GetSymbol() == "N" and atom.IsInRing():
                                    # Check if this N has a methyl group attached
                                    for neighbor in atom.GetNeighbors():
                                        if (
                                            neighbor.GetSymbol() == "C"
                                            and neighbor.GetDegree() == 1
                                        ):
                                            # This might be a methyl group - check if it was added in this reaction
                                            print(
                                                f"Found potential N-CH3 in product at depth {depth}"
                                            )

                                            # Check if this N-CH3 was already present in reactants
                                            n_ch3_in_reactants = False
                                            for reactant in reactants:
                                                if "[N:10]" in reactant and "[CH3]" in reactant:
                                                    n_ch3_in_reactants = True
                                                    break

                                            if not n_ch3_in_reactants:
                                                print(
                                                    f"Late-stage N-alkylation detected at depth {depth} (direct check)"
                                                )
                                                print(f"Reaction SMILES: {rsmi}")
                                                n_alkylation_detected = True
                                                return
                    except Exception as e:
                        print(f"Error in direct N-alkylation check: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_alkylation_detected
