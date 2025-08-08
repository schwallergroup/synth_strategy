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
    This function detects a late-stage functionalization strategy where a significant
    structural modification occurs in the final synthetic steps.
    """
    # Initialize tracking variables
    has_late_stage_functionalization = False

    # List of functional groups commonly used in late-stage functionalization
    late_stage_fg = [
        "Isocyanate",
        "Acyl halide",
        "Sulfonyl halide",
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Azide",
        "Boronic acid",
        "Boronic ester",
        "Nitro group",
        "Diazo",
        "Nitrile",
        "Aldehyde",
        "Ketone",
        "Carboxylic acid",
        "Ester",
        "Anhydride",
        "Halide",
    ]

    # List of reaction types commonly used in late-stage functionalization
    late_stage_reactions = [
        "Suzuki coupling",
        "Buchwald-Hartwig",
        "Sonogashira",
        "Heck",
        "Click chemistry",
        "Acylation of Nitrogen Nucleophiles",
        "Aromatic fluorination",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "N-arylation",
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Minisci",
        "Chan-Lam",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_functionalization

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Focus on final steps (low depth in retrosynthesis)
            try:
                # Extract reactants and products
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a late-stage functionalization reaction
                for reaction_type in late_stage_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected late-stage functionalization reaction: {reaction_type} at depth {depth}"
                        )
                        has_late_stage_functionalization = True
                        return

                # Check for functional group transformations
                for fg in late_stage_fg:
                    # Check if any reactant has the functional group
                    reactants_have_fg = any(checker.check_fg(fg, r) for r in reactants_smiles if r)

                    # Check if the product doesn't have the functional group (transformation occurred)
                    # or if the product has the functional group but reactants don't (addition occurred)
                    product_has_fg = (
                        checker.check_fg(fg, product_smiles) if product_smiles else False
                    )

                    if (reactants_have_fg and not product_has_fg) or (
                        product_has_fg and not reactants_have_fg
                    ):
                        print(f"Detected late-stage functionalization with {fg} at depth {depth}")
                        has_late_stage_functionalization = True
                        return

                # Check for ring formations or modifications
                common_rings = [
                    "benzene",
                    "pyridine",
                    "pyrrole",
                    "furan",
                    "thiophene",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "pyrazole",
                    "triazole",
                    "tetrazole",
                ]

                for ring in common_rings:
                    reactants_have_ring = any(
                        checker.check_ring(ring, r) for r in reactants_smiles if r
                    )
                    product_has_ring = (
                        checker.check_ring(ring, product_smiles) if product_smiles else False
                    )

                    if (reactants_have_ring and not product_has_ring) or (
                        product_has_ring and not reactants_have_ring
                    ):
                        print(f"Detected late-stage ring modification with {ring} at depth {depth}")
                        has_late_stage_functionalization = True
                        return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage functionalization strategy present: {has_late_stage_functionalization}")
    return has_late_stage_functionalization
