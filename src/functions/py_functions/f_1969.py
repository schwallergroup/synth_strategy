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
    This function detects a strategy involving sequential functionalization
    of heteroatoms (N, O, S) through multiple steps.
    """
    heteroatom_functionalization_steps = []
    reaction_types_detected = set()

    # List of reaction types that involve heteroatom functionalization
    heteroatom_functionalization_reactions = [
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Acylation of secondary amines with anhydrides",
        "Alkylation of amines",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "Methylation with MeI_SH",
        "Methylation with DMS",
        "DMS Amine methylation",
        "N-methylation",
        "S-methylation",
        "O-methylation",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "Williamson Ether Synthesis",
        "S-alkylation of thiols",
        "S-alkylation of thiols (ethyl)",
        "S-alkylation of thiols with alcohols",
        "S-alkylation of thiols with alcohols (ethyl)",
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Schotten-Baumann to ester",
        "Esterification of Carboxylic Acids",
        "Urea synthesis via isocyanate and primary amine",
        "Urea synthesis via isocyanate and secondary amine",
        "Urea synthesis via isocyanate and diazo",
        "Urea synthesis via isocyanate and sulfonamide",
    ]

    # List of functional groups containing heteroatoms
    heteroatom_fgs = [
        "Amine",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aniline",
        "Amide",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Ester",
        "Ether",
        "Alcohol",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Phenol",
        "Carboxylic acid",
        "Sulfonamide",
        "Sulfone",
        "Sulfoxide",
        "Thiol",
        "Aromatic thiol",
        "Aliphatic thiol",
        "Urea",
        "Thiourea",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a heteroatom functionalization reaction
                is_functionalization = False
                reaction_type = None

                # Method 1: Check for specific reaction types
                for rxn_type in heteroatom_functionalization_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_functionalization = True
                        reaction_type = rxn_type
                        print(f"Detected reaction type: {rxn_type} at depth {depth}")
                        break

                # Method 2: Check for heteroatom functional group changes
                if not is_functionalization:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check for new heteroatom functional groups in product
                        product_fgs = set()
                        for fg in heteroatom_fgs:
                            if checker.check_fg(fg, product):
                                product_fgs.add(fg)

                        # Check reactants for these functional groups
                        reactant_fgs = set()
                        for reactant in reactants:
                            for fg in heteroatom_fgs:
                                if checker.check_fg(fg, reactant):
                                    reactant_fgs.add(fg)

                        # If there are new functional groups in the product
                        new_fgs = product_fgs - reactant_fgs
                        if new_fgs:
                            is_functionalization = True
                            reaction_type = f"Formation of {', '.join(new_fgs)}"
                            print(
                                f"Detected new functional groups: {new_fgs} at depth {depth}"
                            )

                if is_functionalization:
                    heteroatom_functionalization_steps.append(depth)
                    if reaction_type:
                        reaction_types_detected.add(reaction_type)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If we have at least 3 heteroatom functionalization steps
    if len(heteroatom_functionalization_steps) >= 3:
        print(
            f"Sequential heteroatom functionalization strategy detected at depths: {heteroatom_functionalization_steps}"
        )
        print(f"Reaction types detected: {reaction_types_detected}")
        return True

    print(
        f"Only found {len(heteroatom_functionalization_steps)} heteroatom functionalization steps, need at least 3"
    )
    return False
