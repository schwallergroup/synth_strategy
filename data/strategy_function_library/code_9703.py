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


AMINE_FGS = ["Primary amine", "Secondary amine", "Tertiary amine"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a specific functional group transformation sequence:
    thioether → sulfone → amine substitution.
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

    # Track the transformations in sequence
    transformations = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal findings_json
        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfone to amine substitution
                sulfone_in_reactants = any(checker.check_fg("Sulfone", r) for r in reactants)
                amine_in_reactants = any(any(checker.check_fg(fg, r) for fg in AMINE_FGS) for r in reactants)
                amine_in_product = any(checker.check_fg(fg, product) for fg in AMINE_FGS)
                sulfone_not_in_product = not checker.check_fg("Sulfone", product)

                if sulfone_in_reactants and amine_in_reactants:
                    if amine_in_product and sulfone_not_in_product:
                        print(
                            f"Found potential sulfone to amine substitution at depth {depth}: {rsmi}"
                        )
                        transformations.append(("sulfone_to_amine", depth, rsmi))
                        if "sulfone_to_amine" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("sulfone_to_amine")
                        if "Sulfone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfone")
                        for fg in AMINE_FGS:
                            if checker.check_fg(fg, product) and fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                # Check for thioether oxidation to sulfone
                monosulfide_in_reactants = any(checker.check_fg("Monosulfide", r) for r in reactants)
                sulfone_in_product = checker.check_fg("Sulfone", product)
                sulfone_not_in_reactants = not any(checker.check_fg("Sulfone", r) for r in reactants)

                if monosulfide_in_reactants and sulfone_in_product:
                    if sulfone_not_in_reactants:
                        print(
                            f"Found potential thioether to sulfone oxidation at depth {depth}: {rsmi}"
                        )
                        transformations.append(("thioether_to_sulfone", depth, rsmi))
                        if "thioether_to_sulfone" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("thioether_to_sulfone")
                        if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
                        if "Sulfone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Sulfone")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        elif node["type"] == "mol":
            # Check for thioether in molecule nodes
            if checker.check_fg("Monosulfide", node["smiles"]):
                print(f"Found thioether at depth {depth}: {node['smiles']}")
                transformations.append(("thioether", depth, node["smiles"]))
                if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")

            # Also check for sulfone in molecule nodes (for completeness)
            if checker.check_fg("Sulfone", node["smiles"]):
                print(f"Found sulfone at depth {depth}: {node['smiles']}")
                transformations.append(("sulfone", depth, node["smiles"]))
                if "Sulfone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfone")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth, current_path)
            else: # Assuming 'mol' or other types that should increase depth
                dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Check if we have all required transformations
    has_thioether = any(t[0] == "thioether" for t in transformations)
    has_thioether_to_sulfone = any(t[0] == "thioether_to_sulfone" for t in transformations)
    has_sulfone_to_amine = any(t[0] == "sulfone_to_amine" for t in transformations)

    # Check if transformations are in the correct sequence (by depth)
    in_sequence = False
    if has_thioether and has_thioether_to_sulfone and has_sulfone_to_amine:
        # Get depths for each transformation
        thioether_depths = [t[1] for t in transformations if t[0] == "thioether"]
        sulfone_depths = [t[1] for t in transformations if t[0] == "thioether_to_sulfone"]
        amine_depths = [t[1] for t in transformations if t[0] == "sulfone_to_amine"]

        # Check if there's at least one valid sequence (thioether at higher depth than sulfone, which is at higher depth than amine)
        for t_depth in thioether_depths:
            for s_depth in sulfone_depths:
                for a_depth in amine_depths:
                    if t_depth > s_depth > a_depth:
                        in_sequence = True
                        print(
                            f"Found complete sequence: thioether ({t_depth}) → sulfone ({s_depth}) → amine ({a_depth})"
                        )
                        break
                if in_sequence:
                    break
            if in_sequence:
                break

    # If we have thioether but missing other transformations, check if they might be implicit
    if has_thioether and not (has_thioether_to_sulfone and has_sulfone_to_amine):
        print("Checking for implicit transformations...")
        # Look for molecules containing sulfone at intermediate depths
        sulfone_mols = [t for t in transformations if t[0] == "sulfone"]
        if sulfone_mols and not has_thioether_to_sulfone:
            print("Found sulfone molecules but no explicit oxidation reaction")
            # Add implicit transformation
            transformations.append(
                ("thioether_to_sulfone", min([t[1] for t in sulfone_mols]), "implicit")
            )
            has_thioether_to_sulfone = True

        # Re-check sequence with implicit transformations
        if has_thioether and has_thioether_to_sulfone and has_sulfone_to_amine:
            thioether_depths = [t[1] for t in transformations if t[0] == "thioether"]
            sulfone_depths = [t[1] for t in transformations if t[0] == "thioether_to_sulfone"]
            amine_depths = [t[1] for t in transformations if t[0] == "sulfone_to_amine"]

            for t_depth in thioether_depths:
                for s_depth in sulfone_depths:
                    for a_depth in amine_depths:
                        if t_depth > s_depth > a_depth:
                            in_sequence = True
                            print(
                                f"Found complete sequence with implicit steps: thioether ({t_depth}) → sulfone ({s_depth}) → amine ({a_depth})"
                            )
                            break
                    if in_sequence:
                        break
                if in_sequence:
                    break

    if in_sequence:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    "Monosulfide",
                    "thioether_to_sulfone",
                    "sulfone_to_amine"
                ],
                "description": "A molecule with a Monosulfide (thioether) must appear earlier in the synthesis (higher depth) than the thioether_to_sulfone reaction, which must appear earlier than the sulfone_to_amine reaction."
            }
        })

    print(f"Sequence detection result: {in_sequence}")
    return in_sequence, findings_json