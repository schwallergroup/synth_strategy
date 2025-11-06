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
    Detects a synthesis strategy that involves changing the saturation state
    of nitrogen-containing heterocycles (e.g., pyridine to piperidine).
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

    # Track if we found the required patterns
    pyridine_to_piperidine_conversion = False
    late_stage_reduction = False

    # Track molecules containing pyridine and piperidine at each depth
    pyridine_molecules = {}
    piperidine_molecules = {}

    def dfs_traverse(node, depth=0):
        nonlocal pyridine_to_piperidine_conversion, late_stage_reduction, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Process reaction node
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for hydrogenation reactions
            is_hydrogenation = checker.check_reaction("Hydrogenation (double to single)", rsmi)
            if is_hydrogenation:
                if "Hydrogenation (double to single)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Hydrogenation (double to single)")

            if not is_hydrogenation:
                # Check if any reactant has pyridine and product has piperidine
                product_has_piperidine = checker.check_ring("piperidine", product)
                if product_has_piperidine:
                    if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperidine")

                for reactant in reactants:
                    if checker.check_ring("pyridine", reactant) and product_has_piperidine:
                        if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                        is_hydrogenation = True
                        break

            if is_hydrogenation:
                print(f"Found hydrogenation reaction at depth {depth}: {rsmi}")

                # Check if product contains piperidine and any reactant contains pyridine
                product_has_piperidine = checker.check_ring("piperidine", product)
                if product_has_piperidine:
                    if "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("piperidine")

                reactant_has_pyridine = any(checker.check_ring("pyridine", r) for r in reactants)
                if reactant_has_pyridine:
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                if reactant_has_pyridine and product_has_piperidine:
                    pyridine_to_piperidine_conversion = True
                    if {"type": "sequence", "details": {"description": "A single reaction step must convert a reactant containing a pyridine ring into a product containing a piperidine ring.", "before": "pyridine", "after": "piperidine"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "A single reaction step must convert a reactant containing a pyridine ring into a product containing a piperidine ring.", "before": "pyridine", "after": "piperidine"}})
                    print(f"Found pyridine to piperidine conversion at depth {depth}: {rsmi}")

                    # Check if this is a late-stage reduction (depth 0 or 1)
                    if depth <= 1:
                        late_stage_reduction = True
                        if {"type": "positional", "details": {"description": "The pyridine-to-piperidine reduction event must occur at a late stage in the synthesis.", "target": "pyridine_to_piperidine_reduction", "position": "late_stage", "condition": "depth <= 1"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "The pyridine-to-piperidine reduction event must occur at a late stage in the synthesis.", "target": "pyridine_to_piperidine_reduction", "position": "late_stage", "condition": "depth <= 1"}})
                        print(f"Found late-stage reduction at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases for next node (reaction)
                new_depth = depth + 1
            # If current node is reaction, depth remains the same for next node (chemical)
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # If we didn't find a direct conversion in a reaction, check if the route contains
    # both pyridine and piperidine at different depths
    if not pyridine_to_piperidine_conversion and pyridine_molecules and piperidine_molecules:
        # Check if piperidine appears at an earlier (lower) depth than pyridine
        min_piperidine_depth = min(piperidine_molecules.values())
        for pyridine_depth in pyridine_molecules.values():
            if pyridine_depth > min_piperidine_depth:
                pyridine_to_piperidine_conversion = True
                if {"type": "sequence", "details": {"description": "A single reaction step must convert a reactant containing a pyridine ring into a product containing a piperidine ring.", "before": "pyridine", "after": "piperidine"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "A single reaction step must convert a reactant containing a pyridine ring into a product containing a piperidine ring.", "before": "pyridine", "after": "piperidine"}})
                print(
                    f"Inferred pyridine to piperidine conversion: piperidine at depth {min_piperidine_depth}, pyridine at depth {pyridine_depth}"
                )

                # Check if piperidine is at late stage (depth 0 or 1)
                if min_piperidine_depth <= 1:
                    late_stage_reduction = True
                    if {"type": "positional", "details": {"description": "The pyridine-to-piperidine reduction event must occur at a late stage in the synthesis.", "target": "pyridine_to_piperidine_reduction", "position": "late_stage", "condition": "depth <= 1"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"description": "The pyridine-to-piperidine reduction event must occur at a late stage in the synthesis.", "target": "pyridine_to_piperidine_reduction", "position": "late_stage", "condition": "depth <= 1"}})
                    print(
                        f"Inferred late-stage reduction with piperidine at depth {min_piperidine_depth}"
                    )

    # The strategy is present if we found pyridine to piperidine conversion
    # and it involves a late-stage reduction
    strategy_present = pyridine_to_piperidine_conversion and late_stage_reduction
    print(f"Strategy detection result: {strategy_present}")
    print(f"- Pyridine to piperidine conversion: {pyridine_to_piperidine_conversion}")
    print(f"- Late-stage reduction: {late_stage_reduction}")

    return strategy_present, findings_json
