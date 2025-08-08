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
    Detects a synthetic strategy where the final step involves formation of an amide bond.
    """
    # Initialize tracking variables
    has_amide_formation = False
    late_stage_amide_formation = False

    def dfs_traverse(node, depth=0, is_first_call=True):
        nonlocal has_amide_formation, late_stage_amide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation reactions directly
            is_amide_formation = (
                checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                or checker.check_reaction("Ester with ammonia to amide", rsmi)
                or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
            )

            # If not a known reaction type, check for functional group changes
            if not is_amide_formation:
                # Check for amine in reactants
                reactant_has_amine = any(
                    checker.check_fg("Primary amine", r_smiles)
                    or checker.check_fg("Secondary amine", r_smiles)
                    for r_smiles in reactants_smiles
                )

                # Check for acyl source in reactants
                reactant_has_acyl = any(
                    checker.check_fg("Acyl halide", r_smiles)
                    or checker.check_fg("Carboxylic acid", r_smiles)
                    or checker.check_fg("Ester", r_smiles)
                    or checker.check_fg("Anhydride", r_smiles)
                    for r_smiles in reactants_smiles
                )

                # Check for amide in product
                product_has_amide = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                # If reactants have amine and acyl source, and product has amide, it's likely amide formation
                is_amide_formation = reactant_has_amine and reactant_has_acyl and product_has_amide

            if is_amide_formation:
                print(f"Found amide formation reaction at depth {depth}: {rsmi}")
                has_amide_formation = True

                # In retrosynthesis, late-stage reactions are at depth 0 or 1
                if is_first_call or depth <= 1:
                    late_stage_amide_formation = True
                    print(f"Found late-stage amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, False)

    # Start traversal
    dfs_traverse(route, is_first_call=True)

    # The strategy requires amide formation at a late stage (depth 0 or 1)
    strategy_present = late_stage_amide_formation

    print(f"Late-stage amide formation strategy detected: {strategy_present}")
    return strategy_present
