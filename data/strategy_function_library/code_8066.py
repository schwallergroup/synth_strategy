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


from rdkit import Chem

# Helper functions and classes (e.g., checker) are assumed to be defined elsewhere

RING_FORMING_REACTIONS = [
    "Paal-Knorr pyrrole synthesis",
    "Formation of NOS Heterocycles",
    "Benzothiazole formation",
    "Benzoxazole formation",
    "Benzimidazole formation",
    "Diels-Alder",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
]

COMMON_RINGS_FOR_FORMATION_CHECK = [
    "pyrrole", "pyridine", "furan", "thiophene", "benzene", "indole",
    "pyrazole", "imidazole", "oxazole", "thiazole", "pyrimidine",
    "quinoline", "isoquinoline", "benzothiazole", "benzoxazole", "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step synthetic strategy composed of three key stages: 1) A late-stage Suzuki coupling (depth <= 1). 2) A mid-stage borylation (depth > 1 and <= 3). 3) An early-stage ring formation (depth >= 2). The Suzuki and borylation steps are identified using specific named reaction checkers. The ring formation is identified either by a named reaction from the `RING_FORMING_REACTIONS` list or by the de novo formation of a ring from the `COMMON_RINGS_FOR_FORMATION_CHECK` list. The overall strategy is only confirmed if all three stages are found in the synthesis.
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

    has_suzuki_coupling = False
    has_borylation = False
    has_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling, has_borylation, has_ring_formation, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling (late stage, depth <= 1)
                if depth <= 1:
                    suzuki_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic acids OTf",
                        "Suzuki coupling with boronic esters",
                        "Suzuki coupling with boronic esters OTf"
                    ]
                    for rxn_name in suzuki_reactions:
                        if checker.check_reaction(rxn_name, rsmi):
                            has_suzuki_coupling = True
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            break

                # Check for borylation (mid stage, depth > 1 and <= 3)
                if depth > 1 and depth <= 3:
                    borylation_reactions = [
                        "Preparation of boronic acids",
                        "Preparation of boronic ethers",
                        "Preparation of boronic acids from trifluoroborates"
                    ]
                    for rxn_name in borylation_reactions:
                        if checker.check_reaction(rxn_name, rsmi):
                            has_borylation = True
                            if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                            break

                # Check for ring formation (early stage, depth >= 2)
                if depth >= 2:
                    ring_forming_reaction_found = False
                    for rxn in RING_FORMING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            has_ring_formation = True
                            ring_forming_reaction_found = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break
                    
                    product_rings = [
                        ring for ring in COMMON_RINGS_FOR_FORMATION_CHECK if checker.check_ring(ring, product)
                    ]
                    if product_rings:
                        for ring in product_rings:
                            is_new_ring = True
                            for r in reactants:
                                try:
                                    if Chem.MolFromSmiles(r) and checker.check_ring(ring, r):
                                        is_new_ring = False
                                        break
                                except Exception:
                                    # Handle cases where SMILES might be invalid
                                    pass
                            if is_new_ring:
                                has_ring_formation = True
                                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = has_suzuki_coupling and has_borylation and has_ring_formation

    if has_suzuki_coupling and depth <= 1: # This condition is implicitly handled by the dfs_traverse, but explicitly adding for structural constraint
        # Add positional constraint for late_stage_suzuki_coupling
        if {"type": "positional", "details": {"event_name": "late_stage_suzuki_coupling", "targets": ["Suzuki coupling with boronic acids", "Suzuki coupling with boronic acids OTf", "Suzuki coupling with boronic esters", "Suzuki coupling with boronic esters OTf"], "condition": "depth <= 1"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"event_name": "late_stage_suzuki_coupling", "targets": ["Suzuki coupling with boronic acids", "Suzuki coupling with boronic acids OTf", "Suzuki coupling with boronic esters", "Suzuki coupling with boronic esters OTf"], "condition": "depth <= 1"}})

    if has_borylation and depth > 1 and depth <= 3: # This condition is implicitly handled by the dfs_traverse, but explicitly adding for structural constraint
        # Add positional constraint for mid_stage_borylation
        if {"type": "positional", "details": {"event_name": "mid_stage_borylation", "targets": ["Preparation of boronic acids", "Preparation of boronic ethers", "Preparation of boronic acids from trifluoroborates"], "condition": "1 < depth <= 3"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"event_name": "mid_stage_borylation", "targets": ["Preparation of boronic acids", "Preparation of boronic ethers", "Preparation of boronic acids from trifluoroborates"], "condition": "1 < depth <= 3"}})

    if has_ring_formation and depth >= 2: # This condition is implicitly handled by the dfs_traverse, but explicitly adding for structural constraint
        # Add positional constraint for early_stage_ring_formation
        if {"type": "positional", "details": {"event_name": "early_stage_ring_formation", "targets": ["Paal-Knorr pyrrole synthesis", "Formation of NOS Heterocycles", "Benzothiazole formation", "Benzoxazole formation", "Benzimidazole formation", "Diels-Alder", "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "ring_formation"], "condition": "depth >= 2"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"event_name": "early_stage_ring_formation", "targets": ["Paal-Knorr pyrrole synthesis", "Formation of NOS Heterocycles", "Benzothiazole formation", "Benzoxazole formation", "Benzimidazole formation", "Diels-Alder", "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "ring_formation"], "condition": "depth >= 2"}})

    if result:
        # Add co-occurrence constraint if all three stages are found
        if {"type": "co-occurrence", "details": {"targets": ["late_stage_suzuki_coupling", "mid_stage_borylation", "early_stage_ring_formation"], "description": "The overall strategy requires the presence of all three defined stages: a late-stage Suzuki, a mid-stage borylation, and an early-stage ring formation."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["late_stage_suzuki_coupling", "mid_stage_borylation", "early_stage_ring_formation"], "description": "The overall strategy requires the presence of all three defined stages: a late-stage Suzuki, a mid-stage borylation, and an early-stage ring formation."}})

    return result, findings_json
