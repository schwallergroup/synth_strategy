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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy that involves the formation of a cyclopropane ring and a late-stage amide coupling.
    The amide coupling is identified from a specific list of named reactions. A reaction is considered 'late-stage'
    if it occurs in the final 75% of the synthetic steps.
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

    # Initialize flags for each strategy component
    has_cyclopropane = False
    has_late_amide_coupling = False

    # Track depth for late-stage determination
    max_depth = 0
    amide_coupling_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_cyclopropane, has_late_amide_coupling, max_depth, amide_coupling_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for cyclopropane ring
            if checker.check_ring("cyclopropane", mol_smiles):
                has_cyclopropane = True
                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check for cyclopropanation reaction
            if checker.check_ring("cyclopropane", product) and not any(
                checker.check_ring("cyclopropane", r) for r in reactants
            ):
                has_cyclopropane = True
                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Check for amide coupling reaction from the curated list
            for name in AMIDE_COUPLING_REACTIONS:
                if checker.check_reaction(name, rxn_smiles):
                    # If this is the first amide coupling or it's at a lower depth (later stage)
                    if amide_coupling_depth is None or depth < amide_coupling_depth:
                        amide_coupling_depth = depth
                    if name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                    break # Found one, no need to check other names for this reaction

        # Recursively process children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical)
                new_depth = depth + 1
            # If current node is 'reaction', new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Determine if amide coupling is late-stage (lower depth)
    if amide_coupling_depth is not None:
        # Consider it late-stage if it's in the final 75% of the synthesis depth
        if amide_coupling_depth <= max_depth * 0.75:
            has_late_amide_coupling = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "amide_coupling",
                    "position": "within_final_75_percent_of_steps"
                }
            })

    # Check if both strategy components are present
    combined_strategy = has_cyclopropane and has_late_amide_coupling

    if has_cyclopropane:
        # This constraint is met if has_cyclopropane is True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "cyclopropane",
                "operator": ">=",
                "value": 1
            }
        })

    return combined_strategy, findings_json
