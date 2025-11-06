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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the route involves a late-stage SNAr coupling between a halopyridine
    and a secondary amine (like piperidine).
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

    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected, findings_json

        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_halopyridine = False
                has_secondary_amine = False

                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        if checker.check_ring("pyridine", reactant):
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                            has_halopyridine = True

                    if checker.check_fg("Secondary amine", reactant):
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        has_secondary_amine = True

                if has_halopyridine and has_secondary_amine:
                    # Add co-occurrence constraint
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Aromatic halide on pyridine reactant",
                                "Secondary amine reactant"
                            ],
                            "scope": "reaction_step"
                        }
                    })

                    has_pyridine_in_product = checker.check_ring("pyridine", product)
                    if has_pyridine_in_product:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                    has_aromatic_halide_in_product = checker.check_fg("Aromatic halide", product)
                    if has_aromatic_halide_in_product:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                    has_tertiary_amine_in_product = checker.check_fg("Tertiary amine", product)
                    if has_tertiary_amine_in_product:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

                    if (
                        has_pyridine_in_product
                        and not has_aromatic_halide_in_product
                        and has_tertiary_amine_in_product
                    ):
                        snar_detected = True
                        # Add SNAr reaction and positional constraint
                        findings_json["atomic_checks"]["named_reactions"].append("SNAr")
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "SNAr",
                                "position": "late_stage"
                            }
                        })
                        # Add negation constraint
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "Aromatic halide",
                                "scope": "product_of_SNAr"
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return snar_detected, findings_json
