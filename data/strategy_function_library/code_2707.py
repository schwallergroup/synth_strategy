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


LATE_STAGE_N_INCORPORATION_FGS = [
    "Primary amine",
    "Secondary amine",
    "Tertiary amine",
    "Aniline",
    "Azide",
    "Nitro group",
    "Nitrile",
    "Primary amide",
    "Secondary amide",
    "Tertiary amide",
    "Hydrazine",
    "Hydrazone",
]

N_ARYLATION_REACTIONS = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
    "Ullmann-Goldberg Substitution amine",
    "Buchwald-Hartwig",
    "Ullmann condensation",
    "N-arylation_heterocycles",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where specific nitrogen-containing functional groups are incorporated in the late stages (depth <= 3) of a synthesis. This is identified by checking for the presence of an aromatic halide and a nitrogen nucleophile (from the LATE_STAGE_N_INCORPORATION_FGS list) reacting via a known N-arylation reaction (from the N_ARYLATION_REACTIONS list). The function ensures the nitrogen group is not already present on the aromatic halide reactant.
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

    n_incorporations = 0
    n_incorporation_depths = []

    # Define the full strategy JSON for easy lookup of structural constraints
    full_strategy_json = {
      "function_id": "code_2707",
      "filepath": "../data/merged_good_perf/code_2707.py",
      "description": "Detects a strategy where specific nitrogen-containing functional groups are incorporated in the late stages (depth <= 3) of a synthesis. This is identified by checking for the presence of an aromatic halide and a nitrogen nucleophile reacting via a known N-arylation reaction. The function ensures the nitrogen group is not already present on the aromatic halide reactant.",
      "atomic_checks": {
        "named_reactions": [
          "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
          "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
          "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
          "Goldberg coupling",
          "Goldberg coupling aryl amine-aryl chloride",
          "Goldberg coupling aryl amide-aryl chloride",
          "Ullmann-Goldberg Substitution amine",
          "Buchwald-Hartwig",
          "Ullmann condensation",
          "N-arylation_heterocycles"
        ],
        "ring_systems": [],
        "functional_groups": [
          "Primary amine",
          "Secondary amine",
          "Tertiary amine",
          "Aniline",
          "Azide",
          "Nitro group",
          "Nitrile",
          "Primary amide",
          "Secondary amide",
          "Tertiary amide",
          "Hydrazine",
          "Hydrazone",
          "Aromatic halide"
        ]
      },
      "structural_constraints": [
        {
          "type": "co-occurrence",
          "details": {
            "targets": [
              "any N_ARYLATION_REACTIONS",
              "Aromatic halide",
              "any LATE_STAGE_N_INCORPORATION_FGS"
            ],
            "scope": "reaction_step"
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "N-arylation reaction",
            "position": "late_stage",
            "constraint": "depth <= 3"
          }
        },
        {
          "type": "negation",
          "details": {
            "target": "any LATE_STAGE_N_INCORPORATION_FGS",
            "scope": "Aromatic halide reactant"
          }
        },
        {
          "type": "count",
          "details": {
            "target": "late-stage N-incorporation event",
            "operator": ">=",
            "value": 1
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal n_incorporations, n_incorporation_depths, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                n_nucleophile = False
                n_nucleophile_smiles = None
                halogenated_aromatic = False
                halogenated_reactant = None

                for reactant_smiles in reactants_smiles:
                    for n_fg in LATE_STAGE_N_INCORPORATION_FGS:
                        if checker.check_fg(n_fg, reactant_smiles):
                            n_nucleophile = True
                            n_nucleophile_smiles = reactant_smiles
                            if n_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(n_fg)

                    if checker.check_fg("Aromatic halide", reactant_smiles):
                        halogenated_aromatic = True
                        halogenated_reactant = reactant_smiles
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                is_n_arylation = False
                for rxn_type in N_ARYLATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_n_arylation = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if n_nucleophile and halogenated_aromatic and is_n_arylation and depth <= 3:
                    # Check for negation constraint: N group not already present on aromatic halide reactant
                    negation_satisfied = True
                    if halogenated_reactant:
                        for n_fg in LATE_STAGE_N_INCORPORATION_FGS:
                            if checker.check_fg(n_fg, halogenated_reactant):
                                negation_satisfied = False
                                break

                    if negation_satisfied:
                        if n_nucleophile_smiles != halogenated_reactant:
                            n_incorporations += 1
                            n_incorporation_depths.append(depth)

                            # Record structural constraints if all conditions met for an event
                            # Co-occurrence
                            co_occurrence_constraint = next((c for c in full_strategy_json["structural_constraints"] if c["type"] == "co-occurrence"), None)
                            if co_occurrence_constraint and co_occurrence_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(co_occurrence_constraint))

                            # Positional (depth <= 3)
                            positional_constraint = next((c for c in full_strategy_json["structural_constraints"] if c["type"] == "positional"), None)
                            if positional_constraint and positional_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(positional_constraint))

                            # Negation
                            negation_constraint = next((c for c in full_strategy_json["structural_constraints"] if c["type"] == "negation"), None)
                            if negation_constraint and negation_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(negation_constraint))

            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_present = n_incorporations >= 1

    if strategy_present:
        # Count constraint
        count_constraint = next((c for c in full_strategy_json["structural_constraints"] if c["type"] == "count"), None)
        if count_constraint and count_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(copy.deepcopy(count_constraint))

    return strategy_present, findings_json
