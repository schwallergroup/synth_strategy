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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Carboxylic acid to amide conversion",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the route uses a late-stage amide formation strategy
    (amide formation occurs in the last 2 steps of the synthesis).
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

    has_late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_formation, findings_json

        # Process current node
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is an amide formation reaction using the checker function
            is_amide_formation = False
            for rxn_type in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_amide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            # If not detected by reaction type, check for amide formation by examining product and reactants
            if not is_amide_formation:
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if amide is in product
                has_amide_in_product = False
                for fg_name in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                    if checker.check_fg(fg_name, product):
                        has_amide_in_product = True
                        findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                
                # Check if amide is not in any reactant
                has_amide_in_reactants = False
                for reactant in reactants:
                    for fg_name in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                        if checker.check_fg(fg_name, reactant):
                            has_amide_in_reactants = True
                            # Do not add to findings_json here as it's a negative check
                            break
                    if has_amide_in_reactants:
                        break

                # Check for reactants that could form amides
                has_acid_or_acyl = False
                acid_acyl_fgs = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                for reactant in reactants:
                    for fg_name in acid_acyl_fgs:
                        if checker.check_fg(fg_name, reactant):
                            has_acid_or_acyl = True
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                            break
                    if has_acid_or_acyl:
                        break

                has_amine = False
                amine_fgs = ["Primary amine", "Secondary amine", "Aniline"]
                for reactant in reactants:
                    for fg_name in amine_fgs:
                        if checker.check_fg(fg_name, reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                            break
                    if has_amine:
                        break

                # Amide formation if: amide in product + not in reactants + appropriate reactants present
                is_amide_formation = (
                    has_amide_in_product
                    and not has_amide_in_reactants
                    and has_acid_or_acyl
                    and has_amine
                )

                if is_amide_formation:
                    print(f"Detected amide formation by functional group analysis at depth {depth}")
                    print(f"  Product has amide: {has_amide_in_product}")
                    print(f"  Reactants have amide: {has_amide_in_reactants}")
                    print(f"  Reactants have acid/acyl: {has_acid_or_acyl}")
                    print(f"  Reactants have amine: {has_amine}")
                    # Add a generic 'amide_formation' to named_reactions if detected by FG analysis
                    if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

            # Check if this is a late-stage reaction (depth <= 2 for last 2 steps)
            if is_amide_formation and depth <= 2:
                has_late_amide_formation = True
                print(f"Late-stage amide formation detected at depth {depth}")
                # Add the structural constraint
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "amide_formation",
                        "position": "within_last_three_stages"
                    }
                })

        # Traverse children
        for child in node.get("children", []):
            # The new rule is that depth should only increase when traversing from a chemical node to a reaction node.
            # The depth should remain the same when traversing from a reaction node to a chemical node.
            next_depth = depth
            if node["type"] == "chemical": # Depth increases only when moving from chemical to reaction
                next_depth = depth + 1
            dfs_traverse(child, next_depth)

    # Start traversal from depth 0 (target molecule)
    dfs_traverse(route)
    print(f"Late-stage amide formation: {has_late_amide_formation}")
    
    # Ensure unique entries in lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))
    # Structural constraints are objects, so a simple set won't work. Assuming they are added uniquely by logic.

    return has_late_amide_formation, findings_json
