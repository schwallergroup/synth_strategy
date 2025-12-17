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

HETEROCYCLES_OF_INTEREST = [
    "isoxazole", "furan", "pyrrole", "thiophene", "pyridine", "pyrazole",
    "imidazole", "oxazole", "thiazole", "triazole", "tetrazole", "pyrimidine",
    "pyrazine", "pyridazine", "benzoxazole", "benzothiazole", "benzimidazole",
    "indole", "quinoline", "isoquinoline",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{benzothiazole}", "{benzimidazole_derivatives_aldehyde}",
    "{benzimidazole_derivatives_carboxylic-acid/ester}", "{thiazole}",
    "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}", "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}", "{pyrazole}", "{oxadiazole}",
    "{benzofuran}", "{benzothiophene}", "{indole}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage introduction of specific heterocycles in the final synthetic step. It identifies reactions where the heterocycle is either formed (e.g., via cyclization) or coupled as a pre-existing fragment. The heterocycles and formation reactions checked are defined in `HETEROCYCLES_OF_INTEREST` and `HETEROCYCLE_FORMATION_REACTIONS`.
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
    result = False

    reactions_with_depth = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]
            reactions_with_depth.append((depth, reactants_smiles, product_smiles, rsmi))
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same (new_depth = depth)
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if not reactions_with_depth:
        return False, findings_json

    reactions_with_depth.sort(key=lambda x: x[0])
    final_step = reactions_with_depth[0]
    _, final_step_reactants_smiles, final_step_product_smiles, final_step_rsmi = final_step

    heterocycle_in_final_product = False
    heterocycle_found = None
    for heterocycle in HETEROCYCLES_OF_INTEREST:
        if checker.check_ring(heterocycle, final_step_product_smiles):
            heterocycle_in_final_product = True
            heterocycle_found = heterocycle
            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
            break

    if not heterocycle_in_final_product:
        return False, findings_json

    heterocycle_in_reactants = any(
        checker.check_ring(heterocycle_found, r) for r in final_step_reactants_smiles
    )

    # Case 1: Heterocycle is FORMED in the final step
    if not heterocycle_in_reactants:
        for formation in HETEROCYCLE_FORMATION_REACTIONS:
            if checker.check_reaction(formation, final_step_rsmi):
                result = True
                if formation not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(formation)
                # Add structural constraints for formation
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "heterocycle_introduction_or_formation", "position": "last_stage"}})
                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "specific_heterocycle_formation_reaction"]}})
                return result, findings_json

    # Case 2: Heterocycle is INTRODUCED via coupling in the final step
    if heterocycle_in_reactants:
        reactants_with_heterocycle = [
            r for r in final_step_reactants_smiles if checker.check_ring(heterocycle_found, r)
        ]
        # True if it's a coupling reaction where exactly one reactant brings the heterocycle
        if len(reactants_with_heterocycle) == 1 and len(final_step_reactants_smiles) > 1:
            result = True
            # Add structural constraints for introduction via coupling
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "heterocycle_introduction_or_formation", "position": "last_stage"}})
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_with_target_heterocycle", "operator": "==", "value": 1}})
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "total_reactants_in_step", "operator": ">", "value": 1}})
            return result, findings_json

    return result, findings_json
