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
    This function detects a Friedel-Crafts type acylation pattern in the final step.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and depth <= 1:  # Final or near-final step
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Direct check for Friedel-Crafts acylation reaction
                if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                    found_pattern = True
                    print(f"Found Friedel-Crafts acylation reaction at depth {depth}")
                    return

                # Check for related reactions that might be Friedel-Crafts type acylations
                related_reactions = [
                    "Friedel-Crafts alkylation",
                    "Acylation of olefines by aldehydes",
                    "Friedel-Crafts alkylation with halide",
                ]

                for reaction in related_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        found_pattern = True
                        print(f"Found related reaction {reaction} at depth {depth}")
                        return

                # Fallback to component analysis if the direct checks fail
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if we have an acyl compound and an aromatic compound as reactants
                has_acyl_compound = False
                has_aromatic_compound = False

                for reactant in reactants:
                    # Check for acyl compounds (expanded list)
                    acyl_groups = [
                        "Acyl halide",
                        "Anhydride",
                        "Carboxylic acid",
                        "Ester",
                        "Aldehyde",
                    ]
                    if any(checker.check_fg(fg, reactant) for fg in acyl_groups):
                        has_acyl_compound = True
                        print(f"Found acyl compound at depth {depth}: {reactant}")

                    # Check for aromatic compounds - expanded list
                    aromatic_rings = [
                        "benzene",
                        "naphthalene",
                        "anthracene",
                        "furan",
                        "pyrrole",
                        "thiophene",
                        "pyridine",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "pyrazine",
                        "pyridazine",
                        "triazole",
                        "tetrazole",
                        "carbazole",
                    ]

                    # Try to detect aromatic compounds
                    if any(checker.check_ring(ring, reactant) for ring in aromatic_rings):
                        has_aromatic_compound = True
                        print(f"Found aromatic compound at depth {depth}: {reactant}")
                    elif (
                        Chem.MolFromSmiles(reactant) is not None
                        and Chem.MolFromSmiles(reactant).GetNumAromaticRings() > 0
                    ):
                        has_aromatic_compound = True
                        print(
                            f"Found aromatic compound (RDKit detection) at depth {depth}: {reactant}"
                        )

                # Check if product has a carbonyl group attached to an aromatic ring
                carbonyl_groups = ["Ketone", "Aldehyde"]
                has_carbonyl = any(checker.check_fg(fg, product_part) for fg in carbonyl_groups)
                has_aromatic_product = any(
                    checker.check_ring(ring, product_part) for ring in aromatic_rings
                )

                if has_carbonyl and has_aromatic_product:
                    # If we have both acyl and aromatic reactants, and the product has both carbonyl and aromatic groups
                    if has_acyl_compound and has_aromatic_compound:
                        found_pattern = True
                        print(f"Found Friedel-Crafts type acylation pattern at depth {depth}")

                # Check for specific acylation reactions
                acylation_reactions = [
                    "Acylation of secondary amines with anhydrides",
                    "Acylation of secondary amines",
                    "Acylation of primary amines",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                ]

                for reaction in acylation_reactions:
                    if checker.check_reaction(reaction, rsmi):
                        # For these reactions, verify they involve an aromatic ring
                        if has_aromatic_compound and has_aromatic_product:
                            found_pattern = True
                            print(
                                f"Found acylation reaction {reaction} involving aromatic compounds at depth {depth}"
                            )
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Friedel-Crafts type acylation: {found_pattern}")
    return found_pattern
