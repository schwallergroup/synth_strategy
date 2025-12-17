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


COMPLEX_FRAGMENT_FGS = [
    "Primary amine", "Secondary amine", "Tertiary amine", "Amide", "Ester",
    "Carboxylic acid", "Alcohol", "Ether", "Nitrile", "Nitro group", "Halide",
    "Ketone", "Aldehyde", "Phenol",
]
COMPLEX_FRAGMENT_RINGS = [
    "benzene", "pyridine", "furan", "thiophene", "pyrrole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "piperidine", "morpholine",
    "cyclohexane", "cyclopentane",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis uses a convergent approach where
    complex fragments are combined in late-stage reactions.
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

    late_stage_fragment_count = 0
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_fragment_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            depth = current_depth

            # Focus on late-stage reactions (traversal depth <= 1)
            if depth <= 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "complex_fragment_check",
                        "position": "late_stage (depth <= 1)"
                    }
                })
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            current_ring_checks = []
                            for ring in COMPLEX_FRAGMENT_RINGS:
                                if checker.check_ring(ring, r):
                                    current_ring_checks.append(ring)
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                            ring_count = len(current_ring_checks)

                            current_fg_checks = []
                            for fg in COMPLEX_FRAGMENT_FGS:
                                if checker.check_fg(fg, r):
                                    current_fg_checks.append(fg)
                                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                            fg_count = len(current_fg_checks)

                            is_complex = (
                                (ring_count >= 1 and fg_count >= 1)
                                or ring_count >= 2
                                or fg_count >= 3
                            )

                            if is_complex:
                                late_stage_fragment_count += 1
                    except:
                        # Silently ignore SMILES parsing errors
                        pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    if late_stage_fragment_count >= 2:
        result = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "late_stage_complex_fragment",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
