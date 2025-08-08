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
    Detects if the synthetic route maintains a halogen handle (e.g., aryl bromide)
    through early steps for late-stage functionalization.
    """
    # Track halogen presence through steps
    halogen_at_depth = {}
    halogen_used_at_depth = None
    max_depth = 0

    # Track consecutive depths with halogens in each branch
    branch_paths = []

    # Expanded list of halogen functional groups
    halogen_fgs = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide",
        "Haloalkyne",
    ]

    # Expanded list of reactions that consume halogens
    halogen_consuming_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with sulfonic esters",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_aryl OTf",
        "Stille reaction_vinyl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Negishi coupling",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Heck_non-terminal_vinyl",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_alkenyl halide",
        "Sonogashira acetylene_alkenyl halide",
        "Buchwald-Hartwig",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Goldberg coupling",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Kumada cross-coupling",
        "Hiyama-Denmark Coupling",
        "Aryllithium cross-coupling",
    ]

    def track_consecutive_halogens(node, depth=0, path=None):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if path is None:
            path = []

        new_path = path.copy()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_halogen = any(checker.check_fg(fg, mol_smiles) for fg in halogen_fgs)
            halogen_at_depth[depth] = halogen_at_depth.get(depth, False) or has_halogen

            if has_halogen:
                print(f"Found halogen at depth {depth} in molecule: {mol_smiles}")
                new_path.append(depth)

        # If leaf node or no children, record the path
        if not node.get("children", []):
            if new_path:  # Only record paths that have halogens
                branch_paths.append(new_path)
        else:
            for child in node.get("children", []):
                track_consecutive_halogens(child, depth + 1, new_path)

    def dfs_traverse(node, depth=0):
        nonlocal halogen_used_at_depth, max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for any halogen functional group
            has_halogen = any(checker.check_fg(fg, mol_smiles) for fg in halogen_fgs)
            halogen_at_depth[depth] = halogen_at_depth.get(depth, False) or has_halogen

            if has_halogen:
                print(f"Found halogen at depth {depth} in molecule: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants have halogen but product doesn't (halogen was consumed)
            reactant_has_halogen = any(
                any(checker.check_fg(fg, r) for fg in halogen_fgs) for r in reactants
            )
            product_has_halogen = any(checker.check_fg(fg, product) for fg in halogen_fgs)

            # Check if this is a reaction that typically consumes halogens
            is_coupling_reaction = False
            for reaction_type in halogen_consuming_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    is_coupling_reaction = True
                    break

            if reactant_has_halogen and (not product_has_halogen or is_coupling_reaction):
                if halogen_used_at_depth is None or depth < halogen_used_at_depth:
                    halogen_used_at_depth = depth
                    print(f"Halogen used at depth {depth} in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    track_consecutive_halogens(route)

    # Find the longest consecutive sequence of depths with halogens
    longest_consecutive = 0
    for path in branch_paths:
        path.sort()  # Sort depths in ascending order
        current_consecutive = 1
        max_consecutive_in_path = 1

        for i in range(1, len(path)):
            if path[i] == path[i - 1] + 1:  # Consecutive depths
                current_consecutive += 1
                max_consecutive_in_path = max(max_consecutive_in_path, current_consecutive)
            else:
                current_consecutive = 1

        longest_consecutive = max(longest_consecutive, max_consecutive_in_path)

    print(f"Longest consecutive steps with halogen: {longest_consecutive}")

    # Check if halogen was preserved through multiple steps before being used
    if halogen_used_at_depth is None:
        # Check if halogen is present but never used (might be a valid strategy too)
        halogen_present = any(halogen_at_depth.values())
        if halogen_present:
            print("Halogen present but not used in coupling reactions")
            return longest_consecutive >= 3
        else:
            print("No halogen detected in the route")
            return False

    # Count how many steps the halogen was preserved before being used
    # In retrosynthesis, we need to count from max_depth down to halogen_used_at_depth
    preserved_steps = 0
    for d in range(halogen_used_at_depth + 1, max_depth + 1):
        if d in halogen_at_depth and halogen_at_depth[d]:
            preserved_steps += 1

    print(
        f"Preserved steps: {preserved_steps}, Used at depth: {halogen_used_at_depth}, Max depth: {max_depth}"
    )

    # Return True if halogen was preserved for at least 2 steps before being used
    # and the usage is in a late stage (in retrosynthesis, this means low depth number)
    return preserved_steps >= 2 and halogen_used_at_depth <= 2
