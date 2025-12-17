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


SNAR_REACTION_TYPES = [
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
]

RING_FORMING_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "Intramolecular transesterification/Lactone formation",
    "Intramolecular amination (heterocycle formation)",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the strategy of sequential SNAr reactions on a trifluoromethyl-substituted pyridine
    followed by cyclization.
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

    # Track if we've found each component of the strategy
    cf3_pyridine_nodes = []
    snar_reactions = []
    cyclization_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal cf3_pyridine_nodes, snar_reactions, cyclization_reactions, findings_json

        if node["type"] == "mol":
            # Check for CF3-pyridine in molecule nodes
            mol_smiles = node["smiles"]
            has_cf3 = checker.check_fg("Trifluoro group", mol_smiles)
            has_pyridine = checker.check_ring("pyridine", mol_smiles)

            if has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
            if has_pyridine:
                if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")

            if has_cf3 and has_pyridine:
                cf3_pyridine_nodes.append((depth, mol_smiles))
                # Structural constraint: A molecule containing both a trifluoromethyl group and a pyridine ring
                # This constraint is met if cf3_pyridine_nodes is not empty later.

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for SNAr reactions
                is_snar = False
                for rxn_type in SNAR_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_snar = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_snar:
                    # Verify it's happening on a CF3-pyridine
                    reactants = reactants_smiles.split(".")
                    for r in reactants:
                        if checker.check_fg("Trifluoro group", r) and checker.check_ring(
                            "pyridine", r
                        ):
                            snar_reactions.append((depth, rsmi))
                            # Structural constraint: At least one SNAr reaction must occur on a CF3-pyridine substrate.
                            # This constraint is met if snar_reactions is not empty later.
                            break

                # Check for specific ring-forming reactions
                is_cyclization = False
                for rxn_type in RING_FORMING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_cyclization = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break
                
                if is_cyclization:
                    # Verify it involves a CF3-pyridine
                    if checker.check_ring("pyridine", product_smiles) and checker.check_fg(
                        "Trifluoro group", product_smiles
                    ):
                        cyclization_reactions.append((depth, rsmi))
                        # This cyclization is relevant for the sequence constraint.

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth to analyze the sequence
    cf3_pyridine_nodes.sort(key=lambda x: x[0])
    snar_reactions.sort(key=lambda x: x[0])
    cyclization_reactions.sort(key=lambda x: x[0])

    # Check if we have the complete strategy
    has_cf3_pyridine = len(cf3_pyridine_nodes) > 0
    has_snar = len(snar_reactions) >= 1
    has_cyclization = len(cyclization_reactions) > 0

    # Add structural constraints to findings_json
    if has_cf3_pyridine:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Trifluoro group",
                    "pyridine"
                ],
                "scope": "molecule_in_route",
                "description": "A molecule containing both a trifluoromethyl group and a pyridine ring must be present in the route."
            }
        })
    
    if has_snar:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "SNAr reaction"
                ],
                "scope": "route",
                "description": "At least one SNAr reaction must occur in the route on a trifluoromethyl-pyridine substrate."
            }
        })

    # Check for the sequence - more flexible approach
    correct_sequence = False
    sequence_condition_met = False
    count_condition_met = False

    if has_snar and has_cyclization:
        min_snar_depth = min([r[0] for r in snar_reactions])
        min_cyclization_depth = min([r[0] for r in cyclization_reactions])

        # In retrosynthetic analysis, cyclization should be at a lower depth than at least one SNAr
        if min_cyclization_depth <= min_snar_depth:
            correct_sequence = True
            sequence_condition_met = True

    # If we have multiple SNAr reactions but no detected cyclization, it might still be valid
    if len(snar_reactions) >= 2:
        correct_sequence = True
        count_condition_met = True

    if sequence_condition_met or count_condition_met:
        # Add the alternative_conditions constraint if either sub-condition is met
        alt_conditions_constraint = {
            "type": "alternative_conditions",
            "details": {
                "description": "The route must satisfy one of the following two conditions.",
                "conditions": []
            }
        }
        if sequence_condition_met:
            alt_conditions_constraint["details"]["conditions"].append({
                "type": "sequence",
                "details": {
                    "first_event": "SNAr reaction",
                    "second_event": "cyclization_reaction",
                    "relationship": "before_in_forward_synthesis",
                    "description": "A cyclization reaction must occur at the same synthetic stage or after an SNAr reaction."
                }
            })
        if count_condition_met:
            alt_conditions_constraint["details"]["conditions"].append({
                "type": "count",
                "details": {
                    "target": "SNAr reaction",
                    "operator": ">=",
                    "value": 2,
                    "description": "There must be at least two SNAr reactions in the route."
                }
            })
        findings_json["structural_constraints"].append(alt_conditions_constraint)

    # Check if the complete strategy is present
    strategy_present = (
        has_cf3_pyridine and has_snar and correct_sequence
    )

    return strategy_present, findings_json