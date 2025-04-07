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
    Detects N-oxidation as a late-stage modification (low depth in synthetic tree).
    """
    n_oxidation_at_low_depth = False

    # List of potential oxidation reactions that might include N-oxidation
    oxidation_reactions = [
        "Aromatic hydroxylation",
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation of amide to carboxylic acid",
        "Oxidation of boronic acids",
        "Oxidation of boronic esters",
    ]

    # Nitrogen-containing functional groups that might be involved in N-oxidation
    n_containing_fgs = [
        "Nitro group",
        "Nitroso",
        "Azoxy",
        "Azide",
        "Diazene",
        "Diazo",
        "Nitrate",
        "Nitrile",
        "Isocyanide",
        "Isocyanate",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
    ]

    # Nitrogen-containing rings that might undergo N-oxidation
    n_containing_rings = [
        "pyridine",
        "pyrazole",
        "imidazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "piperidine",
        "piperazine",
        "morpholine",
        "quinoline",
        "isoquinoline",
        "purine",
        "indole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal n_oxidation_at_low_depth

        if node["type"] == "reaction" and depth <= 2:  # Consider depth 0-2 as late-stage
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an oxidation reaction
                is_oxidation = any(checker.check_reaction(rxn, rsmi) for rxn in oxidation_reactions)

                if is_oxidation:
                    print(f"Found potential oxidation reaction at depth {depth}: {rsmi}")

                    # Check if product has nitrogen-containing functional groups
                    product_n_fgs = [fg for fg in n_containing_fgs if checker.check_fg(fg, product)]

                    # Check if product has nitrogen-containing rings
                    product_n_rings = [
                        ring for ring in n_containing_rings if checker.check_ring(ring, product)
                    ]

                    if product_n_fgs or product_n_rings:
                        # Check if there's a change in nitrogen-containing functional groups
                        for reactant in reactants:
                            reactant_n_fgs = [
                                fg for fg in n_containing_fgs if checker.check_fg(fg, reactant)
                            ]
                            if set(product_n_fgs) != set(reactant_n_fgs):
                                print(f"Found N-oxidation (FG change) at depth {depth}: {rsmi}")
                                n_oxidation_at_low_depth = True
                                break

                        # Check if there's a change in nitrogen-containing rings
                        for reactant in reactants:
                            reactant_n_rings = [
                                ring
                                for ring in n_containing_rings
                                if checker.check_ring(ring, reactant)
                            ]
                            if set(product_n_rings) != set(reactant_n_rings):
                                print(f"Found N-oxidation (ring change) at depth {depth}: {rsmi}")
                                n_oxidation_at_low_depth = True
                                break

                # Additional check: Look for specific N-oxidation patterns
                # Check for increase in oxygen count in nitrogen-containing molecules
                if not n_oxidation_at_low_depth:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check if product contains nitrogen
                            has_nitrogen = any(
                                atom.GetAtomicNum() == 7 for atom in product_mol.GetAtoms()
                            )

                            if has_nitrogen:
                                # Count oxygens in product
                                product_o_count = sum(
                                    1 for atom in product_mol.GetAtoms() if atom.GetAtomicNum() == 8
                                )

                                # Check if any reactant has fewer oxygens
                                for reactant in reactants:
                                    reactant_mol = Chem.MolFromSmiles(reactant)
                                    if reactant_mol:
                                        # Check if reactant contains nitrogen
                                        reactant_has_n = any(
                                            atom.GetAtomicNum() == 7
                                            for atom in reactant_mol.GetAtoms()
                                        )

                                        if reactant_has_n:
                                            # Count oxygens in reactant
                                            reactant_o_count = sum(
                                                1
                                                for atom in reactant_mol.GetAtoms()
                                                if atom.GetAtomicNum() == 8
                                            )

                                            # If product has more oxygens than reactant, it might be N-oxidation
                                            if product_o_count > reactant_o_count:
                                                print(
                                                    f"Found potential N-oxidation (O count increase) at depth {depth}: {rsmi}"
                                                )
                                                n_oxidation_at_low_depth = True
                                                break
                    except Exception as e:
                        print(f"Error analyzing molecule: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return n_oxidation_at_low_depth
