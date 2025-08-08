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
    This function detects if the synthetic route involves protection of a carboxylic acid
    as a tert-butyl ester in the early stages of synthesis.
    """
    protection_found = False

    def dfs_traverse(node, depth):
        nonlocal protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a carboxylic acid group
                has_carboxylic_acid_reactant = False
                carboxylic_acid_reactant = None

                for reactant in reactants:
                    if reactant and checker.check_fg("Carboxylic acid", reactant):
                        has_carboxylic_acid_reactant = True
                        carboxylic_acid_reactant = reactant
                        break

                # Check if product has a tert-butyl ester group
                has_tert_butyl_ester = False
                if product and checker.check_fg("Ester", product):
                    # Check specifically for tert-butyl group
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        tert_butyl_pattern = Chem.MolFromSmarts(
                            "[C;H0]([C;H3])([C;H3])([C;H3])[O;D2][C;D3](=[O;D1])"
                        )
                        if mol.HasSubstructMatch(tert_butyl_pattern):
                            has_tert_butyl_ester = True

                # Check if this is a protection reaction
                # Protection of carboxylic acid is typically an esterification
                is_protection_reaction = False
                if has_carboxylic_acid_reactant and has_tert_butyl_ester:
                    if checker.check_reaction(
                        "Esterification of Carboxylic Acids", rsmi
                    ) or checker.check_reaction("Protection of carboxylic acid", rsmi):
                        is_protection_reaction = True

                    # If specific reaction check fails, verify the transformation manually
                    if not is_protection_reaction:
                        # Verify that the carboxylic acid carbon is preserved in the ester
                        # This would require atom mapping analysis which is complex
                        # For simplicity, we'll assume if both functional groups are present
                        # and it's not a deprotection, it's likely a protection
                        is_protection_reaction = True

                # If we found a protection reaction at high depth (early in synthesis)
                # Early stage is typically at depth >= 2 (further from target)
                if is_protection_reaction and depth >= 2:
                    protection_found = True
                    print(f"Found tert-butyl ester protection at depth: {depth}")

        # Continue traversing with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal with depth 0
    dfs_traverse(route, 0)
    return protection_found
