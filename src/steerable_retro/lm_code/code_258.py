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
    Detects a synthetic strategy involving multiple steps that use hydroxylamine
    as a reagent for functional group transformations.
    """
    # Track reactions involving hydroxylamine at different depths
    hydroxylamine_reactions = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for hydroxylamine in reactants string
                hydroxylamine_patterns = [
                    "NH2OH",
                    "ONH2",
                    "[NH2][OH]",
                    "NO",
                    "[NH3+][O-]",
                    "[NH2:",
                    "[OH:",
                ]
                hydroxylamine_present = False

                # Check if any pattern is in the reactants
                if any(pattern in reactants_smiles for pattern in hydroxylamine_patterns):
                    hydroxylamine_present = True
                    print(f"Found hydroxylamine pattern in reaction at depth {depth}")

                # If not found by pattern matching, try using the checker
                if not hydroxylamine_present:
                    try:
                        hydroxylamine_present = checker.check_fg("Hydroxylamine", reactants_smiles)
                        if hydroxylamine_present:
                            print(f"Found hydroxylamine using checker at depth {depth}")
                    except:
                        # If Hydroxylamine is not in the fg_dict, check individual reactants
                        reactants = reactants_smiles.split(".")
                        for reactant in reactants:
                            # Check for patterns in individual reactants
                            if any(pattern in reactant for pattern in hydroxylamine_patterns):
                                hydroxylamine_present = True
                                print(
                                    f"Found hydroxylamine pattern in reactant at depth {depth}: {reactant}"
                                )
                                break

                            # Check for NH2 and OH in the same molecule
                            if ("NH2" in reactant and "OH" in reactant) or (
                                "[NH2" in reactant and "[OH" in reactant
                            ):
                                hydroxylamine_present = True
                                print(f"Found NH2 and OH in reactant at depth {depth}: {reactant}")
                                break

                            try:
                                if checker.check_fg("Hydroxylamine", reactant):
                                    hydroxylamine_present = True
                                    print(
                                        f"Found hydroxylamine using checker in reactant at depth {depth}: {reactant}"
                                    )
                                    break
                            except:
                                pass

                # If hydroxylamine is present, check for known hydroxylamine-mediated transformations
                if hydroxylamine_present:
                    print(f"Hydroxylamine present in reaction at depth {depth}")
                    # Check for specific hydroxylamine-mediated reactions
                    if checker.check_reaction(
                        "N-hydroxyimidamide from nitrile and hydroxylamine", rsmi
                    ):
                        print(f"Found N-hydroxyimidamide formation at depth {depth}")
                        hydroxylamine_reactions.add(depth)
                    elif checker.check_reaction("Amidoxime from nitrile and hydroxylamine", rsmi):
                        print(f"Found amidoxime formation at depth {depth}")
                        hydroxylamine_reactions.add(depth)
                    elif checker.check_reaction(
                        "1,2,4-oxadiazol-5(2H)-one synthesis from nitrile, hydrogen carbonate, and hydroxylamine",
                        rsmi,
                    ):
                        print(f"Found 1,2,4-oxadiazol-5(2H)-one synthesis at depth {depth}")
                        hydroxylamine_reactions.add(depth)
                    # Check for carbonyl to oxime transformation
                    elif checker.check_fg("Aldehyde", reactants_smiles) and checker.check_fg(
                        "Oxime", products_smiles
                    ):
                        print(f"Found aldehyde to oxime transformation at depth {depth}")
                        hydroxylamine_reactions.add(depth)
                    elif checker.check_fg("Ketone", reactants_smiles) and checker.check_fg(
                        "Oxime", products_smiles
                    ):
                        print(f"Found ketone to oxime transformation at depth {depth}")
                        hydroxylamine_reactions.add(depth)
                    else:
                        # General case: hydroxylamine is present but specific reaction not identified
                        print(f"Found hydroxylamine usage in reaction at depth {depth}")
                        hydroxylamine_reactions.add(depth)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Total hydroxylamine-mediated transformations found: {len(hydroxylamine_reactions)}")
    print(f"Hydroxylamine reactions found at depths: {sorted(list(hydroxylamine_reactions))}")

    # Return True if hydroxylamine is used in multiple steps (different depths)
    return len(hydroxylamine_reactions) >= 2
