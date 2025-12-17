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


# Refactored for Enumeration
AROMATIC_RINGS_OF_INTEREST = [
    "benzene", "pyridine", "pyrrole", "furan", "thiophene", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "indole", "quinoline",
    "isoquinoline", "naphthalene", "benzoxazole", "benzothiazole", "benzimidazole"
]
NON_AROMATIC_PRECURSORS_OF_INTEREST = ["cyclopentane", "cyclohexane", "cyclopropane"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step sequence where a non-aromatic ring is first formed, followed by a subsequent aromatization step. The initial ring formation is checked against the list `NON_AROMATIC_PRECURSORS_OF_INTEREST`. The aromatization is confirmed by detecting a dehydrogenation reaction, an increase in aromatic atoms, or the formation of a corresponding aromatic ring from the `AROMATIC_RINGS_OF_INTEREST` list.
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

    # Original strategy JSON for structural constraints lookup
    original_strategy_json = {
      "function_id": "code_3824",
      "filepath": "../data/merged_good_perf/code_3824.py",
      "description": "Detects a two-step sequence where a non-aromatic ring is first formed, followed by a subsequent aromatization step. The initial ring formation is checked against a list of non-aromatic precursors. The aromatization is confirmed by detecting a dehydrogenation reaction, an increase in aromatic atoms, or the formation of a corresponding aromatic ring.",
      "atomic_checks": {
        "named_reactions": [
          "ring_formation",
          "Quinone formation",
          "Dehydrogenation"
        ],
        "ring_systems": [
          "benzene",
          "pyridine",
          "pyrrole",
          "furan",
          "thiophene",
          "imidazole",
          "oxazole",
          "thiazole",
          "pyrimidine",
          "indole",
          "quinoline",
          "isoquinoline",
          "naphthalene",
          "benzoxazole",
          "benzothiazole",
          "benzimidazole",
          "cyclopentane",
          "cyclohexane",
          "cyclopropane"
        ],
        "functional_groups": []
      },
      "structural_constraints": [
        {
          "type": "sequence",
          "details": {
            "description": "A non-aromatic ring is formed, which is subsequently aromatized via a Dehydrogenation reaction.",
            "ordered_events": [
              "ring_formation",
              "Dehydrogenation"
            ]
          }
        },
        {
          "type": "sequence",
          "details": {
            "description": "A non-aromatic ring is formed, which is subsequently aromatized via a Quinone formation reaction.",
            "ordered_events": [
              "ring_formation",
              "Quinone formation"
            ]
          }
        },
        {
          "type": "sequence",
          "details": {
            "description": "A non-aromatic ring is formed, followed by a subsequent reaction that forms a new aromatic ring.",
            "ordered_events": [
              "ring_formation",
              "ring_formation"
            ]
          }
        }
      ]
    }

    sequence_found = False
    reactions_by_depth = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactions_by_depth.append((depth, rsmi))
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    reactions_by_depth.sort(key=lambda x: x[0], reverse=True)

    formed_rings_by_reaction = {}

    # Mapping of non-aromatic to potential aromatic counterparts
    aromatic_counterparts = {
        "cyclopentane": ["furan", "pyrrole", "thiophene"],
        "cyclohexane": ["benzene", "pyridine"],
        "cyclopropane": [],
    }

    # Check for ring formation
    for i, (current_depth, current_rsmi) in enumerate(reactions_by_depth):
        try:
            current_reactants = current_rsmi.split(">")[0].split(".")
            current_product = current_rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(current_product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in current_reactants if r]

            if not product_mol or None in reactant_mols:
                continue

            product_rings = set()
            for ring_name in AROMATIC_RINGS_OF_INTEREST + NON_AROMATIC_PRECURSORS_OF_INTEREST:
                if checker.check_ring(ring_name, current_product):
                    product_rings.add(ring_name)
                    if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring_name)

            reactant_rings = set()
            for reactant in current_reactants:
                for ring_name in AROMATIC_RINGS_OF_INTEREST + NON_AROMATIC_PRECURSORS_OF_INTEREST:
                    if checker.check_ring(ring_name, reactant):
                        reactant_rings.add(ring_name)
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)

            formed_rings = product_rings - reactant_rings

            if formed_rings:
                formed_rings_by_reaction[i] = formed_rings
                # Record 'ring_formation' if a new ring is formed
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        except Exception:
            continue

    # Check for aromatization of previously formed rings
    for ring_formation_idx, formed_rings in formed_rings_by_reaction.items():
        # Check subsequent reactions for aromatization
        for j in range(ring_formation_idx + 1, len(reactions_by_depth)):
            next_depth, next_rsmi = reactions_by_depth[j]

            try:
                next_reactants = next_rsmi.split(">")[0].split(".")
                next_product = next_rsmi.split(">")[-1]
                aromatization_detected = False

                # Check for specific, chemically sound aromatization reactions
                oxidation_reactions = ["Quinone formation", "Dehydrogenation"]
                for rxn_name in oxidation_reactions:
                    if checker.check_reaction(rxn_name, next_rsmi):
                        aromatization_detected = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        
                        # Record structural constraint for specific aromatization reactions
                        if rxn_name == "Dehydrogenation":
                            constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("details", {}).get("ordered_events") == ["ring_formation", "Dehydrogenation"]), None)
                            if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)
                        elif rxn_name == "Quinone formation":
                            constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("details", {}).get("ordered_events") == ["ring_formation", "Quinone formation"]), None)
                            if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)
                        break
                if aromatization_detected:
                    sequence_found = True
                    break

                # Check for increase in aromaticity by atom count
                next_product_mol = Chem.MolFromSmiles(next_product)
                next_reactant_mols = [Chem.MolFromSmiles(r) for r in next_reactants if r]

                if next_product_mol and None not in next_reactant_mols:
                    product_aromatic_count = sum(1 for atom in next_product_mol.GetAtoms() if atom.GetIsAromatic())
                    reactant_aromatic_count = sum(sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) for mol in next_reactant_mols)

                    if product_aromatic_count > reactant_aromatic_count:
                        aromatization_detected = True

                    # Check if any of the formed non-aromatic rings now have aromatic counterparts
                    if not aromatization_detected:
                        for ring_name in formed_rings:
                            if ring_name in NON_AROMATIC_PRECURSORS_OF_INTEREST:
                                for aromatic_ring in aromatic_counterparts.get(ring_name, []):
                                    if checker.check_ring(aromatic_ring, next_product) and not any(checker.check_ring(aromatic_ring, r) for r in next_reactants):
                                        aromatization_detected = True
                                        # Record the new aromatic ring formed
                                        if aromatic_ring not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append(aromatic_ring)
                                        # Record 'ring_formation' for the new aromatic ring
                                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                        
                                        # Record structural constraint for new aromatic ring formation
                                        constraint_obj = next((c for c in original_strategy_json["structural_constraints"] if c.get("details", {}).get("ordered_events") == ["ring_formation", "ring_formation"]), None)
                                        if constraint_obj and constraint_obj not in findings_json["structural_constraints"]:
                                            findings_json["structural_constraints"].append(constraint_obj)
                                        break
                            if aromatization_detected:
                                break

                if aromatization_detected:
                    sequence_found = True
                    break

            except Exception:
                continue

        if sequence_found:
            break

    return sequence_found, findings_json