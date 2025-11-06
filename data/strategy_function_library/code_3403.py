from typing import Tuple, Dict, List
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


HETEROCYCLES_OF_INTEREST = [
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "triazole",
    "tetrazole",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "furan",
    "thiophene",
    "morpholine",
    "piperidine",
    "piperazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if specific heterocyclic rings, defined in HETEROCYCLES_OF_INTEREST, that are present in the target molecule are also preserved in late-stage intermediates. Late-stage is defined as intermediates within the first 2 synthetic steps from the target molecule.
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

    target_mol_smiles = route["smiles"]

    target_rings = []
    for ring in HETEROCYCLES_OF_INTEREST:
        if checker.check_ring(ring, target_mol_smiles):
            target_rings.append(ring)
            findings_json["atomic_checks"]["ring_systems"].append(ring)

    if not target_rings:
        # Add structural constraint for 'any_monitored_heterocycle_in_target' not found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_monitored_heterocycle_in_target",
                "operator": ">=",
                "value": 1
            }
        })
        return False, findings_json
    else:
        # Add structural constraint for 'any_monitored_heterocycle_in_target' found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_monitored_heterocycle_in_target",
                "operator": ">=",
                "value": 1
            }
        })

    core_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal core_preserved, findings_json

        if not core_preserved:
            return

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Only check late-stage intermediates (depth 1 and 2).
            # depth=0 is the target, which we don't need to re-check.
            if 0 < depth < 3 and not node.get("in_stock", False):
                all_rings_preserved_in_intermediate = True
                for ring in target_rings:
                    if not checker.check_ring(ring, mol_smiles):
                        core_preserved = False
                        all_rings_preserved_in_intermediate = False
                        # If a ring is lost, record the negation constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "loss_of_target_heterocycle",
                                "context": "intermediate_at_depth_1_or_2"
                            }
                        })
                        return
                # If all rings are preserved in this intermediate, and we haven't added the negation constraint yet
                # This part is tricky as the constraint is about 'loss', so we only add it if loss occurs.
                # If core_preserved remains True after the loop, it means no loss occurred for this intermediate.

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    
    # If core_preserved is still True at the end, it means no loss occurred.
    # The 'negation' constraint is only added if 'loss_of_target_heterocycle' happens.
    # So, if core_preserved is True, we don't add the negation constraint.

    return core_preserved, findings_json
