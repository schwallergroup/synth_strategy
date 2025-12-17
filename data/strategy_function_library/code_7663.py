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


HETEROCYCLIC_RINGS_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane",
    "dioxathiolane", "thiazolidine", "oxazolidine", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole", "benzoxazole",
    "benzothiazole", "benzimidazole", "pteridin", "phenothiazine",
    "phenoxazine", "dibenzofuran", "dibenzothiophene", "xanthene",
    "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine",
    "porphyrin", "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step sequence where a heterocyclic ring is opened in an earlier synthetic step and a heterocyclic ring is closed in a later step. The specific rings considered are defined in the `HETEROCYCLIC_RINGS_OF_INTEREST` list. This identifies strategies where a ring is temporarily opened for modification before being re-formed.
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

    # In retrosynthesis:
    # - ring_formations: rings formed in retrosynthesis (rings broken in forward synthesis)
    # - ring_breakings: rings broken in retrosynthesis (rings formed in forward synthesis)
    ring_formations = []  # (depth, ring_type) -> Forward breaking
    ring_breakings = []  # (depth, ring_type) -> Forward formation

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result
        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    reactant_rings = {}
                    for reactant in reactants:
                        for ring_type in HETEROCYCLIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring_type, reactant):
                                if ring_type not in reactant_rings:
                                    reactant_rings[ring_type] = 0
                                reactant_rings[ring_type] += 1
                                if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_type)

                    product_rings = {}
                    for ring_type in HETEROCYCLIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring_type, product):
                            if ring_type not in product_rings:
                                product_rings[ring_type] = 0
                            product_rings[ring_type] += 1
                            if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_type)

                    # Ring formation in forward direction
                    for ring_type, count in product_rings.items():
                        react_count = reactant_rings.get(ring_type, 0)
                        if count > react_count:
                            ring_breakings.append((depth, ring_type))
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                    # Ring breaking in forward direction
                    for ring_type, count in reactant_rings.items():
                        prod_count = product_rings.get(ring_type, 0)
                        if count > prod_count:
                            ring_formations.append((depth, ring_type))
                            if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                except Exception:
                    pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if ring_formations and ring_breakings:
        # A ring opening (forward breaking) must occur at a greater depth (earlier step)
        # than a ring closing (forward formation, smaller depth, later step).
        for formation_depth, _ in ring_formations:
            for breaking_depth, _ in ring_breakings:
                if formation_depth > breaking_depth:
                    result = True
                    # Add structural constraints if the condition is met
                    co_occurrence_constraint = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "ring_destruction",
                                "ring_formation"
                            ]
                        }
                    }
                    if co_occurrence_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(co_occurrence_constraint)

                    sequence_constraint = {
                        "type": "sequence",
                        "details": {
                            "before": "ring_destruction",
                            "after": "ring_formation"
                        }
                    }
                    if sequence_constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(sequence_constraint)

                    # Since we found one, we can break and return
                    return result, findings_json

    return result, findings_json
