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


THIOETHER_FORMATION_REACTIONS = [
    "thioether_nucl_sub",
    "S-alkylation of thiols",
    "S-alkylation of thiols with alcohols",
    "S-alkylation of thiols (ethyl)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (depth <= 2) thioether formation involving at least one fragment with 8 or more heavy atoms. This strategy is identified by confirming the new formation of a 'Monosulfide' group and checking if the reaction is one of the following types: 'thioether_nucl_sub', 'S-alkylation of thiols', 'S-alkylation of thiols with alcohols', or 'S-alkylation of thiols (ethyl)'.
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

    # Original strategy JSON for structural constraints lookup
    strategy_json_constraints = [
        {
            "type": "positional",
            "details": {
                "target": "thioether_formation",
                "position": "late_stage (depth <= 2)"
            }
        },
        {
            "type": "negation",
            "details": {
                "target": "Monosulfide",
                "scope": "reactants"
            }
        },
        {
            "type": "count",
            "details": {
                "target": "reactant_with_heavy_atoms_gte_8",
                "operator": ">=",
                "value": 1
            }
        }
    ]

    late_stage_thioether_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_thioether_coupling_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth from ID - default to late stage (0) if not specified
            # The depth parameter in dfs_traverse now correctly reflects the current node's depth
            # based on the new traversal rules.

            # Consider depths 0, 1, and 2 as late-stage
            is_late_stage = False
            if depth <= 2:
                is_late_stage = True
                # Record positional constraint if met
                findings_json["structural_constraints"].append(strategy_json_constraints[0])

            # Check if this is a thioether coupling reaction
            is_thioether_reaction = False
            for r_name in THIOETHER_FORMATION_REACTIONS:
                if checker.check_reaction(r_name, rsmi):
                    is_thioether_reaction = True
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)

            # Check if product contains a thioether (monosulfide)
            has_thioether_product = checker.check_fg("Monosulfide", product)
            if has_thioether_product:
                if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")

            # Check if thioether is newly formed (not present in reactants)
            thioether_in_reactants = any(checker.check_fg("Monosulfide", r) for r in reactants)
            if not thioether_in_reactants:
                # Record negation constraint if met
                findings_json["structural_constraints"].append(strategy_json_constraints[1])

            # Check if coupling connects complex fragments
            complex_fragments = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if (
                    mol and mol.GetNumHeavyAtoms() >= 8
                ):  # Define "complex" as having 8+ heavy atoms
                    complex_fragments += 1
            connects_complex_fragments = complex_fragments >= 1
            if connects_complex_fragments:
                # Record count constraint if met
                findings_json["structural_constraints"].append(strategy_json_constraints[2])

            # If all conditions are met, we've found a late-stage thioether coupling
            if (
                is_late_stage
                and is_thioether_reaction
                and has_thioether_product
                and not thioether_in_reactants
                and connects_complex_fragments
            ):
                late_stage_thioether_coupling_found = True

        for child in node.get("children", []):
            # New depth calculation logic:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route, 0)

    # Remove duplicate structural constraints if any were added multiple times
    unique_structural_constraints = []
    seen_constraints = set()
    for constraint in findings_json["structural_constraints"]:
        constraint_str = str(constraint) # Convert dict to string for set comparison
        if constraint_str not in seen_constraints:
            unique_structural_constraints.append(constraint)
            seen_constraints.add(constraint_str)
    findings_json["structural_constraints"] = unique_structural_constraints

    return late_stage_thioether_coupling_found, findings_json
