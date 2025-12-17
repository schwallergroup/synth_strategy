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
    "Ketone",
    "Aldehyde",
    "Ester",
    "Amide",
    "Primary amine",
    "Secondary amine",
    "Tertiary amine",
    "Ether",
    "Alcohol",
    "Carboxylic acid",
    "Nitrile",
    "Nitro group",
    "Halide",
    "Phenol",
    "Aromatic halide",
    "Sulfonamide",
    "Urea",
    "Thiourea",
]

COMPLEX_FRAGMENT_RINGS = [
    "benzene",
    "pyridine",
    "pyrrole",
    "furan",
    "thiophene",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "piperidine",
    "morpholine",
    "cyclohexane",
    "cyclopentane",
    "naphthalene",
    "indole",
    "quinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis pattern in the final step. A step is considered convergent if at least two reactant fragments are deemed 'complex'. A fragment is defined as complex if the sum of its unique ring types (from `COMPLEX_FRAGMENT_RINGS`) and unique functional group types (from `COMPLEX_FRAGMENT_FGS`) is two or more.
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

    print(f"Starting analysis of route with target: {route.get('smiles', 'Unknown')}")
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json

        print(
            f"Traversing node at depth {depth}, type: {node['type']}, smiles: {node.get('smiles', 'N/A')}"
        )

        if node["type"] == "reaction" and depth == 1:
            # Structural constraint: positional - last_stage
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "convergent_synthesis_check",
                    "position": "last_stage"
                }
            })
            try:
                rsmi = node.get("metadata", {}).get("rsmi")
                if rsmi:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Final reaction: {rsmi}")
                    print(f"Number of reactants: {len(reactants)}")

                    if len(reactants) >= 2:
                        # Structural constraint: count - reactants_in_last_stage
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "reactants_in_last_stage",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        complex_fragments = 0

                        for i, reactant in enumerate(reactants):
                            if not reactant:
                                continue

                            print(f"\nAnalyzing reactant {i+1}: {reactant}")

                            ring_count = 0
                            found_rings_for_reactant = []
                            for ring in COMPLEX_FRAGMENT_RINGS:
                                if checker.check_ring(ring, reactant):
                                    ring_count += 1
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    found_rings_for_reactant.append(ring)
                                    print(f"  Found ring: {ring}")

                            fg_count = 0
                            found_fgs_for_reactant = []
                            for fg in COMPLEX_FRAGMENT_FGS:
                                if checker.check_fg(fg, reactant):
                                    fg_count += 1
                                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                                    found_fgs_for_reactant.append(fg)
                                    print(f"  Found functional group: {fg}")

                            print(
                                f"  Ring count: {ring_count}, Functional group count: {fg_count}"
                            )

                            if (ring_count + fg_count) >= 2:
                                # Structural constraint: count - sum_of_unique_rings_and_fgs_per_reactant
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "sum_of_unique_rings_and_fgs_per_reactant",
                                        "operator": ">=",
                                        "value": 2
                                    }
                                })
                                complex_fragments += 1
                                print(f"  Reactant {i+1} is a complex fragment")

                        print(f"\nTotal complex fragments: {complex_fragments}")
                        if complex_fragments >= 2:
                            # Structural constraint: count - complex_reactants_in_last_stage
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "complex_reactants_in_last_stage",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
                            is_convergent = True
                            print(
                                "Detected convergent synthesis with complex fragments in final step"
                            )
            except Exception as e:
                print(f"Error analyzing final reaction: {e}")

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth if node['type'] == 'reaction' else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return is_convergent, findings_json
