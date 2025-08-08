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
    This function detects if the synthetic route involves amide formation in the final step.
    """
    amide_formation_at_depth_zero = False

    def dfs_traverse(node):
        nonlocal amide_formation_at_depth_zero

        if node["type"] == "reaction" and node.get("metadata", {}).get("depth", 0) == 0:
            rsmi = node.get("metadata", {}).get("rsmi")
            if rsmi:
                try:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check if this is an amide formation reaction using the checker
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann to ester",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Carboxylic acid to amide conversion",
                    ]

                    # Check if any of the amide formation reactions match
                    is_amide_formation = any(
                        checker.check_reaction(rxn_name, rsmi)
                        for rxn_name in amide_formation_reactions
                    )

                    # If reaction checker didn't identify it, do a more detailed check
                    if not is_amide_formation:
                        # Check for functional groups in reactants and products
                        acyl_sources = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                        amine_types = ["Primary amine", "Secondary amine", "Aniline"]
                        amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]

                        has_acyl_source = any(
                            any(checker.check_fg(fg, r) for fg in acyl_sources)
                            for r in reactants_smiles
                        )
                        has_amine = any(
                            any(checker.check_fg(amine, r) for amine in amine_types)
                            for r in reactants_smiles
                        )

                        product_has_amide = any(
                            checker.check_fg(amide, product_smiles) for amide in amide_types
                        )

                        # Check if we have the right reactants and product for amide formation
                        if product_has_amide and has_acyl_source and has_amine:
                            # Verify that the amide is formed in this reaction (not pre-existing)
                            reactants_have_amide = any(
                                any(checker.check_fg(amide, r) for amide in amide_types)
                                for r in reactants_smiles
                            )

                            if not reactants_have_amide:
                                is_amide_formation = True

                    if is_amide_formation:
                        amide_formation_at_depth_zero = True

                except Exception as e:
                    pass  # Silently handle errors

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return amide_formation_at_depth_zero
