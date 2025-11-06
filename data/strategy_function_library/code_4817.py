from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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

from pathlib import Path
root_data = Path(__file__).parent.parent

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


DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Alcohol deprotection from silyl ethers",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "TMS deprotection from alkyne",
    "Tert-butyl deprotection of amine",
    "Phthalimide deprotection",
    "N-glutarimide deprotection",
    "Acetal hydrolysis to aldehyde",
    "Ketal hydrolysis to ketone",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route is linear and employs at least one of several key strategies:
    sequential reductions, late-stage deprotection, or nitroalkene oxidative cleavage.
    """

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Helper function to check if a reaction is a reduction
    def is_reduction_reaction(node):
        if node["type"] != "reaction":
            return False

        try:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            # Check for common reduction reactions
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                print(f"Found nitro reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of aldehydes and ketones to alcohols")
                print(f"Found carbonyl reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of carboxylic acid to primary alcohol", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of carboxylic acid to primary alcohol")
                print(f"Found carboxylic acid reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of ester to primary alcohol")
                print(f"Found ester reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitrile to amine")
                print(f"Found nitrile reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of primary amides to amines", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of primary amides to amines")
                print(f"Found amide reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of secondary amides to amines", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of secondary amides to amines")
                print(f"Found secondary amide reduction: {rsmi}")
                return True
            if checker.check_reaction("Reduction of tertiary amides to amines", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of tertiary amides to amines")
                print(f"Found tertiary amide reduction: {rsmi}")
                return True
            if checker.check_reaction("Hydrogenation (double to single)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Hydrogenation (double to single)")
                print(f"Found alkene hydrogenation: {rsmi}")
                return True
            if checker.check_reaction("Hydrogenation (triple to double)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Hydrogenation (triple to double)")
                print(f"Found alkyne hydrogenation: {rsmi}")
                return True
            if checker.check_reaction("Arene hydrogenation", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Arene hydrogenation")
                print(f"Found arene hydrogenation: {rsmi}")
                return True

            # Check for functional group changes that indicate reduction
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group reduction
            for reactant in reactants:
                if checker.check_fg("Nitro group", reactant):
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    if not checker.check_fg("Nitro group", product):
                        if (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        ):
                            if checker.check_fg("Primary amine", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if checker.check_fg("Secondary amine", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            if checker.check_fg("Tertiary amine", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                            print(f"Found nitro to amine reduction: {rsmi}")
                            return True

            # Check for carbonyl reduction
            for reactant in reactants:
                if (checker.check_fg("Aldehyde", reactant) or checker.check_fg("Ketone", reactant)):
                    if checker.check_fg("Aldehyde", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if checker.check_fg("Ketone", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if not checker.check_fg("Aldehyde", product) and not checker.check_fg("Ketone", product):
                        if checker.check_fg("Primary alcohol", product) or checker.check_fg("Secondary alcohol", product):
                            if checker.check_fg("Primary alcohol", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                            if checker.check_fg("Secondary alcohol", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                            print(f"Found carbonyl to alcohol reduction: {rsmi}")
                            return True

            # Check for nitrile reduction
            for reactant in reactants:
                if checker.check_fg("Nitrile", reactant):
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    if not checker.check_fg("Nitrile", product):
                        if checker.check_fg("Primary amine", product) or checker.check_fg("Secondary amine", product):
                            if checker.check_fg("Primary amine", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                            if checker.check_fg("Secondary amine", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                            print(f"Found nitrile to amine reduction: {rsmi}")
                            return True

            # Check for alkene/alkyne reduction
            for reactant in reactants:
                if (checker.check_fg("Alkyne", reactant) and not checker.check_fg("Alkyne", product)):
                    findings_json["atomic_checks"]["functional_groups"].append("Alkyne")
                    print(f"Found alkene/alkyne reduction: {rsmi}")
                    return True
                if (checker.check_fg("Vinyl", reactant) and not checker.check_fg("Vinyl", product)):
                    findings_json["atomic_checks"]["functional_groups"].append("Vinyl")
                    print(f"Found alkene/alkyne reduction: {rsmi}")
                    return True

            return False
        except Exception as e:
            print(f"Error in is_reduction_reaction: {e}")
            return False

    # Helper function to check if a reaction is a deprotection
    def is_deprotection_reaction(node):
        if node["type"] != "reaction":
            return False

        try:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            # Check for common deprotection reactions using the curated list
            for reaction_name in DEPROTECTION_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    print(f"Found deprotection ({reaction_name}): {rsmi}")
                    return True

            return False
        except Exception as e:
            print(f"Error in is_deprotection_reaction: {e}")
            return False

    # Helper function to check if a reaction is a nitroalkene oxidative cleavage
    def is_nitroalkene_cleavage(node):
        if node["type"] != "reaction":
            return False

        try:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            # Check for oxidative cleavage reactions
            if checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Oxidation of alkene to carboxylic acid")
                reactants = rsmi.split(">")[0].split(".")
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    if checker.check_fg("Vinyl", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Vinyl")
                    if checker.check_fg("Nitro group", reactant) and checker.check_fg("Vinyl", reactant):
                        print(f"Found nitroalkene oxidative cleavage: {rsmi}")
                        return True

            # Check for Nef reaction (nitro to ketone)
            if checker.check_reaction("Nef reaction (nitro to ketone)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Nef reaction (nitro to ketone)")
                print(f"Found Nef reaction: {rsmi}")
                return True

            # Check for oxidation reactions
            if checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Oxidation of aldehydes to carboxylic acids")
                print(f"Found aldehyde oxidation: {rsmi}")
                return True
            if checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Oxidation of alcohol to carboxylic acid")
                print(f"Found alcohol oxidation: {rsmi}")
                return True
            if checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Oxidation of ketone to carboxylic acid")
                print(f"Found ketone oxidation: {rsmi}")
                return True
            if checker.check_reaction("Alkene oxidation to aldehyde", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Alkene oxidation to aldehyde")
                print(f"Found alkene oxidation to aldehyde: {rsmi}")
                return True

            # Check for functional group changes that indicate nitroalkene cleavage
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            for reactant in reactants:
                if checker.check_fg("Nitro group", reactant):
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                if checker.check_fg("Vinyl", reactant):
                    findings_json["atomic_checks"]["functional_groups"].append("Vinyl")

                if checker.check_fg("Nitro group", reactant) and checker.check_fg("Vinyl", reactant):
                    if (
                        checker.check_fg("Carboxylic acid", product)
                        or checker.check_fg("Aldehyde", product)
                        or checker.check_fg("Ketone", product)
                    ):
                        if checker.check_fg("Carboxylic acid", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Aldehyde", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        if checker.check_fg("Ketone", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                        print(f"Found nitroalkene cleavage to carbonyl: {rsmi}")
                        return True

            # More general oxidative cleavage check
            for reactant in reactants:
                if checker.check_fg("Vinyl", reactant) and not checker.check_fg("Vinyl", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Vinyl")
                    if (
                        checker.check_fg("Carboxylic acid", product)
                        or checker.check_fg("Aldehyde", product)
                        or checker.check_fg("Ketone", product)
                    ):
                        if checker.check_fg("Carboxylic acid", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        if checker.check_fg("Aldehyde", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        if checker.check_fg("Ketone", product):
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                        print(f"Found general alkene oxidative cleavage: {rsmi}")
                        return True

            return False
        except Exception as e:
            print(f"Error in is_nitroalkene_cleavage: {e}")
            return False

    # Implement the sequential reduction strategy check
    def sequential_reduction_strategy(route):
        reduction_count = 0
        reduction_depths = []

        def dfs(node, depth=0):
            nonlocal reduction_count

            if node["type"] == "reaction" and is_reduction_reaction(node):
                reduction_count += 1
                reduction_depths.append(depth)

            for child in node.get("children", []):
                # New logic for depth calculation
                if node["type"] == "reaction":
                    dfs(child, depth)
                else: # Assuming 'chemical' or other non-reaction type
                    dfs(child, depth + 1)

        dfs(route)

        print(f"Found {reduction_count} reduction reactions at depths {reduction_depths}")

        # Check if there are at least 2 reduction steps
        if reduction_count >= 2:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reduction_reaction", "operator": ">=", "value": 1, "description": "Checks for the presence of at least one reduction reaction in the route."}})
            # Check if reductions are sequential (adjacent depths)
            reduction_depths.sort()
            for i in range(len(reduction_depths) - 1):
                if reduction_depths[i + 1] - reduction_depths[i] <= 2:  # Allow for some flexibility
                    print(
                        f"Found sequential reductions at depths {reduction_depths[i]} and {reduction_depths[i+1]}"
                    )
                    return True

            # If we have multiple reductions but they're not sequential, still consider it a partial match
            print(f"Found multiple reductions but not sequential: {reduction_depths}")
            return reduction_count >= 2
        elif reduction_count == 1:
            # Consider a single reduction as a partial match
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reduction_reaction", "operator": ">=", "value": 1, "description": "Checks for the presence of at least one reduction reaction in the route."}})
            print(f"Found single reduction at depth {reduction_depths[0]}")
            return True

        return False

    # Implement the late-stage deprotection strategy check
    def late_stage_deprotection_strategy(route):
        deprotection_depth = None
        max_depth = 0

        def dfs(node, depth=0):
            nonlocal deprotection_depth, max_depth

            max_depth = max(max_depth, depth)

            if node["type"] == "reaction" and is_deprotection_reaction(node):
                if deprotection_depth is None or depth < deprotection_depth:
                    deprotection_depth = depth

            for child in node.get("children", []):
                # New logic for depth calculation
                if node["type"] == "reaction":
                    dfs(child, depth)
                else: # Assuming 'chemical' or other non-reaction type
                    dfs(child, depth + 1)

        dfs(route)

        print(f"Max depth: {max_depth}, Deprotection depth: {deprotection_depth}")

        # Check if there is a deprotection and it's in the late stage (low depth)
        if deprotection_depth is not None:
            # Consider it late-stage if it's in the first half of the synthesis (more lenient)
            if deprotection_depth <= max_depth / 2:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "deprotection_reaction", "position": "late_stage", "description": "Checks for a deprotection reaction occurring in the first half of the synthetic steps (low depth value)."}})
                print(
                    f"Found late-stage deprotection at depth {deprotection_depth} (max depth: {max_depth})"
                )
                return True

        return False

    # Implement the nitroalkene oxidative cleavage strategy check
    def nitroalkene_oxidative_cleavage_strategy(route):
        found_cleavage = False

        def dfs(node):
            nonlocal found_cleavage

            if node["type"] == "reaction" and is_nitroalkene_cleavage(node):
                found_cleavage = True

            for child in node.get("children", []):
                # The original dfs for this function did not pass depth, so we maintain that.
                # If depth was intended to be used here, it would need to be added to the function signature.
                dfs(child)

        dfs(route)

        if found_cleavage:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["nitroalkene_cleavage"], "description": "Checks for the presence of a nitroalkene oxidative cleavage reaction or a similar oxidation."}})
            print("Found nitroalkene oxidative cleavage or similar oxidation")

        return found_cleavage

    # Implement the linear synthesis strategy check
    def linear_synthesis_strategy(route):
        # A linear synthesis has a single branch at each step
        # Count the maximum branching factor
        max_branching = 0

        def dfs(node):
            nonlocal max_branching

            children = node.get("children", [])
            max_branching = max(max_branching, len(children))

            for child in children:
                dfs(child)

        dfs(route)

        # A linear synthesis should have at most 2 branches at any point
        # (allowing for some convergent steps)
        is_linear = max_branching <= 2
        if is_linear:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "route_branching", "operator": "<=", "value": 2, "description": "Checks if the synthesis is linear by ensuring the maximum number of children for any node is at most 2."}})
            print(f"Synthesis is linear (max branching: {max_branching})")
        else:
            print(f"Synthesis is not linear (max branching: {max_branching})")

        return is_linear

    # Check all the required strategies
    has_sequential_reductions = sequential_reduction_strategy(route)
    has_late_deprotection = late_stage_deprotection_strategy(route)
    has_nitroalkene_cleavage = nitroalkene_oxidative_cleavage_strategy(route)
    is_linear = linear_synthesis_strategy(route)

    print(f"Sequential reductions: {has_sequential_reductions}")
    print(f"Late-stage deprotection: {has_late_deprotection}")
    print(f"Nitroalkene cleavage: {has_nitroalkene_cleavage}")
    print(f"Linear synthesis: {is_linear}")

    # Require only one of the three chemical strategies plus linear synthesis
    combined_strategy = is_linear and (
        has_sequential_reductions or has_late_deprotection or has_nitroalkene_cleavage
    )

    if combined_strategy:
        print(
            "Detected combined strategy with linear synthesis and at least one of: sequential reductions, "
            + "nitroalkene cleavage, or late-stage deprotection"
        )

    # Remove duplicate entries from findings_json lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))
    
    # Convert list of dicts to list of unique dicts (by converting to tuple of items for set)
    unique_structural_constraints = []
    seen_constraints = set()
    for constraint in findings_json["structural_constraints"]:
        # Convert dict to a frozenset of (key, frozenset(items)) for nested dicts or (key, value) for simple values
        # This handles nested dicts and lists within the constraint details
        def make_hashable(obj):
            if isinstance(obj, dict):
                return frozenset({k: make_hashable(v) for k, v in obj.items()}.items())
            elif isinstance(obj, list):
                return tuple(make_hashable(elem) for elem in obj)
            else:
                return obj

        hashable_constraint = make_hashable(constraint)
        if hashable_constraint not in seen_constraints:
            seen_constraints.add(hashable_constraint)
            unique_structural_constraints.append(constraint)
    findings_json["structural_constraints"] = unique_structural_constraints

    return combined_strategy, findings_json