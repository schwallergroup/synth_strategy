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
    "furan", "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole",
    "thiazole", "pyrimidine", "pyrazine", "thiophene", "indole", "quinoline",
    "isoquinoline", "benzimidazole", "benzoxazole", "benzothiazole", "triazole",
    "tetrazole", "oxadiazole", "thiadiazole", "isoxazole", "isothiazole",
]

FORMYLATION_REACTIONS_OF_INTEREST = [
    "Bouveault aldehyde synthesis",
    "Duff reaction",
    "Reimer-Tiemann reaction",
    "Vilsmeier-Haack reaction",
    "Gattermann reaction",
    "Gattermann-Koch reaction"
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects late-stage formylation on a heterocyclic scaffold.
    Late-stage means the formylation occurs in the final synthetic step (depth 1).
    """
    print("Starting late_stage_formylation_strategy analysis")
    formylation_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Original strategy JSON for structural constraints reference
    original_strategy_json = {
      "function_id": "code_394",
      "filepath": "../data/merged_good_perf/code_394.py",
      "description": "This function detects late-stage formylation on a heterocyclic scaffold. Late-stage means the formylation occurs in the final synthetic step (depth 1).",
      "atomic_checks": {
        "named_reactions": [
          "Bouveault aldehyde synthesis",
          "Duff reaction",
          "Reimer-Tiemann reaction",
          "Vilsmeier-Haack reaction",
          "Gattermann reaction",
          "Gattermann-Koch reaction"
        ],
        "ring_systems": [
          "furan",
          "pyrrole",
          "pyridine",
          "pyrazole",
          "imidazole",
          "oxazole",
          "thiazole",
          "pyrimidine",
          "pyrazine",
          "thiophene",
          "indole",
          "quinoline",
          "isoquinoline",
          "benzimidazole",
          "benzoxazole",
          "benzothiazole",
          "triazole",
          "tetrazole",
          "oxadiazole",
          "thiadiazole",
          "isoxazole",
          "isothiazole"
        ],
        "functional_groups": [
          "Aldehyde"
        ]
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "formylation_on_heterocycle",
            "position": "last_stage"
          }
        },
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "any_of_FORMYLATION_REACTIONS_OF_INTEREST",
              "any_of_HETEROCYCLES_OF_INTEREST",
              "Aldehyde"
            ],
            "scope": "same_reaction_step_product"
          }
        },
        {
          "type": "negation",
          "details": {
            "target": "Aldehyde",
            "scope": "reactants_of_formylation_step"
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal formylation_detected, findings_json

        if node.get("children", []):
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # Only increase depth if current node is chemical
                    new_depth = depth + 1
                dfs_traverse(child, new_depth)

        if node["type"] == "reaction" and depth == 1:  # Final step
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # An aldehyde must be formed in the product.
                    aldehyde_in_product = checker.check_fg("Aldehyde", product)
                    if not aldehyde_in_product:
                        return
                    findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")

                    # The product must contain one of the specified heterocycles.
                    has_heterocycle = False
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, product):
                            has_heterocycle = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                    if not has_heterocycle:
                        return

                    # The aldehyde group must be newly formed, not present in reactants.
                    has_aldehyde_in_reactants = any(checker.check_fg("Aldehyde", reactant) for reactant in reactants)

                    # The reaction must be a known, named formylation reaction.
                    is_formylation_reaction = False
                    for rxn in FORMYLATION_REACTIONS_OF_INTEREST:
                        if checker.check_reaction(rxn, rsmi):
                            is_formylation_reaction = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                    if not is_formylation_reaction:
                        return

                    # Final logic: Aldehyde appears on a heterocycle via a known formylation reaction.
                    if not has_aldehyde_in_reactants and is_formylation_reaction:
                        print(f"Late-stage formylation detected on heterocyclic scaffold via {rsmi}")
                        formylation_detected = True

                        # Record structural constraints
                        # Positional constraint: last_stage
                        findings_json["structural_constraints"].append(original_strategy_json["structural_constraints"][0])

                        # Co-occurrence constraint
                        findings_json["structural_constraints"].append(original_strategy_json["structural_constraints"][1])

                        # Negation constraint
                        findings_json["structural_constraints"].append(original_strategy_json["structural_constraints"][2])

                except Exception as e:
                    print(f"Error processing SMILES in formylation detection: {e}")

    dfs_traverse(route)
    print(f"Final result: formylation_detected = {formylation_detected}")
    return formylation_detected, findings_json