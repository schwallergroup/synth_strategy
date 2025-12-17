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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


BOC_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
]

BOC_DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if Boc protection is maintained through multiple steps.
    A proper Boc protection strategy involves:
    1. Adding Boc protection early in the synthesis
    2. Maintaining the Boc group through intermediate steps
    3. Removing the Boc group in a late stage

    Note: In retrosynthetic analysis, higher depth = earlier stage, lower depth = later stage
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

    # Track Boc-protected molecules and protection/deprotection reactions
    protected_molecules = []
    protection_reactions = []
    deprotection_reactions = []

    # Track depth of each reaction for determining early vs late stage
    reaction_depths = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check if molecule has Boc group
            if checker.check_fg("Boc", node["smiles"]):
                protected_molecules.append((depth, node["smiles"]))
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")
                print(f"Found Boc-protected molecule at depth {depth}: {node['smiles'][:30]}...")

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Store reaction depth for later analysis
                reaction_depths[rsmi] = depth

                # Check if this is a Boc protection reaction
                for name in BOC_PROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        protection_reactions.append((depth, rsmi))
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        print(f"Found Boc protection reaction at depth {depth}: {rsmi[:50]}...")
                        break

                # Check if this is a Boc deprotection reaction
                for name in BOC_DEPROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        deprotection_reactions.append((depth, rsmi))
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        print(f"Found Boc deprotection reaction at depth {depth}: {rsmi[:50]}...")
                        break

                # Alternative check: see if Boc group appears or disappears
                has_boc_in_reactants = any(checker.check_fg("Boc", r) for r in reactants)
                has_boc_in_product = checker.check_fg("Boc", product)
                if has_boc_in_reactants or has_boc_in_product:
                    if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")

                # In retrosynthetic analysis, if product has Boc but reactants don't,
                # it's a deprotection reaction when going forward
                if has_boc_in_reactants and not has_boc_in_product:
                    if not any(rsmi == r[1] for r in deprotection_reactions):
                        deprotection_reactions.append((depth, rsmi))
                        print(f"Found implicit Boc deprotection at depth {depth}: {rsmi[:50]}...")

                # In retrosynthetic analysis, if reactants have Boc but product doesn't,
                # it's a protection reaction when going forward
                if has_boc_in_product and not has_boc_in_reactants:
                    if not any(rsmi == r[1] for r in protection_reactions):
                        protection_reactions.append((depth, rsmi))
                        print(f"Found implicit Boc protection at depth {depth}: {rsmi[:50]}...")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Analyze the protection strategy
    print(f"Found {len(protection_reactions)} Boc protection reactions")
    print(f"Found {len(deprotection_reactions)} Boc deprotection reactions")
    print(f"Found {len(protected_molecules)} Boc-protected molecules")

    result = False

    # No protection or deprotection found
    if not protection_reactions and not deprotection_reactions:
        print("No Boc protection/deprotection reactions found")
        return result, findings_json

    # If we found deprotection but no protection, try to infer protection
    if deprotection_reactions and not protection_reactions:
        print("Found deprotection but no explicit protection, checking for protected molecules")
        # If we have protected molecules at deeper levels, infer protection happened
        deep_protected = [
            mol
            for depth, mol in protected_molecules
            if depth > max(d for d, _ in deprotection_reactions)
        ]
        if deep_protected:
            # Infer protection at the deepest level where we see a Boc group
            deepest_protected = max(protected_molecules, key=lambda x: x[0])
            protection_reactions.append(deepest_protected)
            print(f"Inferred Boc protection at depth {deepest_protected[0]}")

    # Still no complete strategy
    if not protection_reactions or not deprotection_reactions:
        print("No complete Boc protection strategy found")
        return result, findings_json

    # If both protection and deprotection are found, add the co-occurrence constraint
    if protection_reactions and deprotection_reactions:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Boc protection",
                    "Boc deprotection"
                ],
                "description": "The route must contain at least one Boc protection reaction and one Boc deprotection reaction. These can be identified by name or by the appearance/disappearance of a Boc group across a reaction step."
            }
        })

    # In retrosynthetic analysis:
    # - Early stage = higher depth
    # - Late stage = lower depth
    early_threshold = max_depth * 2 // 3  # Early stage is in the deeper part of the tree
    late_threshold = max_depth // 3  # Late stage is in the shallower part of the tree

    earliest_protection = max(protection_reactions, key=lambda x: x[0])[
        0
    ]  # Highest depth = earliest
    latest_deprotection = min(deprotection_reactions, key=lambda x: x[0])[
        0
    ]  # Lowest depth = latest

    print(
        f"Earliest protection at depth {earliest_protection}, latest deprotection at depth {latest_deprotection}"
    )
    print(f"Early threshold: {early_threshold}, late threshold: {late_threshold}")

    # Check if there are protected intermediates between protection and deprotection
    # In retrosynthetic analysis, deprotection (lower depth) happens before protection (higher depth)
    intermediate_steps = False
    for depth, mol in protected_molecules:
        if latest_deprotection < depth < earliest_protection:
            intermediate_steps = True
            print(f"Found protected intermediate at depth {depth}")
            break

    # A good Boc strategy has:
    # - Early protection (high depth in retrosynthesis)
    # - Late deprotection (low depth in retrosynthesis)
    # - Protected intermediates between them
    is_good_strategy = (
        earliest_protection >= early_threshold
        and latest_deprotection <= late_threshold
        and intermediate_steps
    )

    if earliest_protection >= early_threshold:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc protection",
                "position": "early_stage",
                "description": "The first Boc protection event must occur in an early stage of the synthesis (retrosynthetically, at a high depth)."
            }
        })

    if latest_deprotection <= late_threshold:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc deprotection",
                "position": "late_stage",
                "description": "The final Boc deprotection event must occur in a late stage of the synthesis (retrosynthetically, at a low depth)."
            }
        })

    if intermediate_steps:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "event_sequence": [
                    "Boc protection",
                    "Boc",
                    "Boc deprotection"
                ],
                "description": "A molecule containing a Boc functional group must exist as an intermediate in the steps between the first protection and the final deprotection."
            }
        })

    result = is_good_strategy
    print(f"Is good Boc protection strategy: {is_good_strategy}")
    return result, findings_json