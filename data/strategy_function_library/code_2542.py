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
    Identifies a synthetic strategy where an early-stage 'Aromatic bromination' is followed by one or more specific late-stage transformations, including 'Acylation of primary amines', 'O-methylation', or 'Benzoxazole formation'.
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

    # Track key transformations and their depths
    bromination_depth = None
    acetylation_depth = None
    o_methylation_depth = None
    ring_formation_depth = None

    result = False # Initialize the main boolean result

    print("Starting analysis of synthesis route...")

    def dfs_traverse(node, current_depth=0):
        nonlocal bromination_depth, acetylation_depth, o_methylation_depth, ring_formation_depth, findings_json

        if node["type"] == "reaction":
            # Get metadata
            metadata = node.get("metadata", {})
            depth = metadata.get("depth", current_depth)
            rsmi = metadata.get("rsmi", "")

            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for bromination
                if checker.check_reaction("Aromatic bromination", rsmi):
                    print(f"Bromination detected at depth: {depth}")
                    bromination_depth = depth
                    if "Aromatic bromination" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Aromatic bromination")

                # Check for acetylation (NH2 → NHAc)
                if checker.check_reaction("Acylation of primary amines", rsmi):
                    print(f"Acetylation detected at depth: {depth}")
                    acetylation_depth = depth
                    if "Acylation of primary amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Acylation of primary amines")

                # Check for O-methylation (OH → OCH3)
                if checker.check_reaction("O-methylation", rsmi):
                    print(f"O-Methylation detected at depth: {depth}")
                    o_methylation_depth = depth
                    if "O-methylation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("O-methylation")

                # Check for benzoxazole formation
                benzoxazole_formation_detected = False
                if checker.check_reaction("Benzoxazole formation from aldehyde", rsmi):
                    print(f"Benzoxazole formation from aldehyde detected at depth: {depth}")
                    ring_formation_depth = depth
                    benzoxazole_formation_detected = True
                    if "Benzoxazole formation from aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Benzoxazole formation from aldehyde")

                if (checker.check_ring("benzoxazole", product_smiles)
                    and not any(checker.check_ring("benzoxazole", r) for r in reactants_smiles)
                ):
                    print(f"Benzoxazole ring formation detected at depth: {depth}")
                    if not benzoxazole_formation_detected: # Avoid double counting if already caught by reaction name
                        ring_formation_depth = depth
                    if "benzoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("benzoxazole")
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else:
                dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Transformation depths - Bromination: {bromination_depth}, Acetylation: {acetylation_depth}, O-Methylation: {o_methylation_depth}, Ring Formation: {ring_formation_depth}"
    )

    # Check if bromination is detected and occurs early (higher depth)
    if bromination_depth is not None:
        # Check if bromination occurs before other transformations
        # Note: Higher depth means earlier in the synthesis (retrosynthetic perspective)
        early_bromination = True

        if acetylation_depth is not None and bromination_depth <= acetylation_depth:
            early_bromination = False

        if o_methylation_depth is not None and bromination_depth <= o_methylation_depth:
            early_bromination = False

        if ring_formation_depth is not None and bromination_depth <= ring_formation_depth:
            early_bromination = False

        # We need at least one other transformation after bromination
        has_subsequent_transformations = (
            acetylation_depth is not None
            or o_methylation_depth is not None
            or ring_formation_depth is not None
        )

        if early_bromination and has_subsequent_transformations:
            print("Linear synthesis with early bromination strategy detected")
            result = True
            # Record structural constraints
            # Constraint: At least one of the subsequent transformations occurred
            if has_subsequent_transformations:
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": [
                            "Acylation of primary amines",
                            "O-methylation",
                            "Benzoxazole formation"
                        ],
                        "operator": ">=",
                        "value": 1
                    }
                })
            # Sequence constraints
            if acetylation_depth is not None and bromination_depth > acetylation_depth:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "Aromatic bromination",
                        "after": "Acylation of primary amines"
                    }
                })
            if o_methylation_depth is not None and bromination_depth > o_methylation_depth:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "Aromatic bromination",
                        "after": "O-methylation"
                    }
                })
            if ring_formation_depth is not None and bromination_depth > ring_formation_depth:
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "Aromatic bromination",
                        "after": "Benzoxazole formation"
                    }
                })

        else:
            print("Bromination detected but not in the expected pattern")
    else:
        print("No bromination detected in the synthesis route")

    return result, findings_json
