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
    This function detects if the synthesis involves a late-stage esterification
    (carboxylic acid → ester in the final steps).

    In retrosynthetic analysis, this means we're looking for ester → acid transformations
    in the late stages of the synthesis route.
    """
    esterification_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_detected

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check late-stage reactions (depth 0, 1, or 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, we're looking for ester in reactants and acid in product
                has_ester = any(checker.check_fg("Ester", r) for r in reactants if r)
                has_acid = product and checker.check_fg("Carboxylic acid", product)

                # Check if this is an esterification or hydrolysis reaction
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )
                is_hydrolysis = checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )
                is_ester_saponification = checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)

                print(f"  Has ester in reactants: {has_ester}")
                print(f"  Has acid in product: {has_acid}")
                print(f"  Is esterification reaction: {is_esterification}")
                print(f"  Is hydrolysis reaction: {is_hydrolysis}")
                print(f"  Is ester saponification: {is_ester_saponification}")

                # In retrosynthesis, a late-stage esterification would appear as a hydrolysis
                # or as an esterification reaction with ester in reactants and acid in product
                if has_ester and has_acid:
                    # This is a potential esterification in retrosynthesis
                    # Check if the reaction is a known esterification/hydrolysis type
                    if is_esterification or is_hydrolysis or is_ester_saponification:
                        esterification_detected = True
                        print(
                            f"Late-stage esterification strategy detected at depth {depth} (known reaction type)"
                        )
                    else:
                        # If not a known reaction type, check if the transformation involves
                        # converting an ester to a carboxylic acid (in retrosynthesis)
                        # by looking at the atom mapping
                        try:
                            # Find a reactant with an ester
                            ester_reactant = next(
                                (r for r in reactants if r and checker.check_fg("Ester", r)), None
                            )
                            if ester_reactant and product:
                                # This is a heuristic approach - if we have an ester converting to an acid
                                # in retrosynthesis, it's likely an esterification strategy
                                esterification_detected = True
                                print(
                                    f"Late-stage esterification strategy detected at depth {depth} (heuristic)"
                                )
                        except Exception as e:
                            print(f"Error in atom mapping analysis: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage esterification strategy detected: {esterification_detected}")
    return esterification_detected
