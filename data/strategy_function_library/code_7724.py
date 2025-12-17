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


# Refactoring for Enumeration: Isolate the list of heterocycles
HETEROCYCLES_OF_INTEREST = [
    "pyrazole", "pyrrole", "furan", "thiophene", "pyridine", "imidazole",
    "oxazole", "thiazole", "triazole", "tetrazole", "pyrimidine", "pyrazine",
    "pyridazine", "indole", "benzimidazole", "morpholine", "piperidine",
    "piperazine", "pyrrolidine", "oxirane", "aziridine", "thiirane",
    "oxetane", "azetidine", "thietane",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis with late-stage amide coupling of fragments
    where one fragment contains a cyclopropane ring and another contains a heterocycle.
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

    amide_coupling_reactions = []  # Store (depth, reactants) tuples

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_reactions, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                is_amide_coupling = False
                amide_reaction_names = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]
                for reaction_name in amide_reaction_names:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amide_coupling = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if is_amide_coupling:
                    reactants = [r for r in rsmi.split(">")[0].split(".") if r]
                    amide_coupling_reactions.append((depth, reactants))

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = False

    if not amide_coupling_reactions:
        return result, findings_json

    amide_coupling_reactions.sort(key=lambda x: x[0])
    latest_depth, latest_reactants = amide_coupling_reactions[0]

    cyclopropane_reactants = []
    heterocycle_reactants = []

    for reactant in latest_reactants:
        if checker.check_ring("cyclopropane", reactant):
            cyclopropane_reactants.append(reactant)
            if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")

        for heterocycle in HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(heterocycle, reactant):
                heterocycle_reactants.append(reactant)
                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                break

    different_fragments = False
    for cp_reactant in cyclopropane_reactants:
        for het_reactant in heterocycle_reactants:
            if cp_reactant != het_reactant:
                different_fragments = True
                break
        if different_fragments:
            break

    # A reaction is at depth=1 in the final step of the synthesis.
    is_late_stage = latest_depth == 1
    is_convergent = len(latest_reactants) >= 2

    if is_late_stage:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Amide Formation",
                "position": "last_stage"
            }
        })
    
    if is_convergent:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reactants",
                "operator": ">=",
                "value": 2,
                "scope": "last_stage_amide_formation"
            }
        })

    if different_fragments:
        # This constraint is only added if both cyclopropane and heterocycle reactants are found
        # AND they are from different fragments.
        if cyclopropane_reactants and heterocycle_reactants:
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "reactant_with_cyclopropane",
                        "reactant_with_heterocycle"
                    ],
                    "scope": "last_stage_amide_formation",
                    "condition": "different_reactants"
                }
            })

    result = is_convergent and different_fragments and is_late_stage

    return result, findings_json