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
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step strategy involving: 1) silyl ether protection of an alcohol, 2) silyl ether deprotection, 3) reduction of a nitro group to an amine, and 4) a late-stage amide coupling.
    The amide coupling is identified by checking for a specific set of reaction types defined in AMIDE_COUPLING_REACTIONS.
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

    # Initialize flags for each strategy
    has_silyl_protection = False
    has_silyl_deprotection = False
    has_nitro_reduction = False
    has_late_amide_coupling = False

    # Track depth for late-stage determination
    max_depth = 0
    amide_coupling_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_silyl_protection, has_silyl_deprotection, has_nitro_reduction
        nonlocal has_late_amide_coupling, max_depth, amide_coupling_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        # Process reaction nodes
        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for silyl protection
            if not has_silyl_protection and checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                has_silyl_protection = True
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")

            # Check for silyl deprotection
            if not has_silyl_deprotection and checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi):
                has_silyl_deprotection = True
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers")

            # Check for nitro reduction
            if not has_nitro_reduction and checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                has_nitro_reduction = True
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

            # Check for amide coupling
            if not has_late_amide_coupling:
                for reaction_type in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        has_late_amide_coupling = True
                        amide_coupling_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

        # Determine next_depth based on current node type
        current_node_type = node.get("type")
        if current_node_type == "reaction":
            next_depth = depth  # Depth remains the same when going from reaction to chemical
        else:
            next_depth = depth + 1  # Depth increases when going from chemical to reaction

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if amide coupling is late-stage (lower depth)
    is_late_stage = False
    if amide_coupling_depth is not None:
        # Consider it late-stage if it's in the first half of the synthesis depth
        is_late_stage = amide_coupling_depth <= max_depth // 2

    # Return True only if all required strategies are present
    combined_strategy = (
        has_silyl_protection
        and has_silyl_deprotection
        and has_nitro_reduction
        and has_late_amide_coupling
        and is_late_stage
    )

    # Add structural constraints if the overall strategy is met
    if combined_strategy:
        # Co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Alcohol protection with silyl ethers",
                    "Alcohol deprotection from silyl ethers",
                    "Reduction of nitro groups to amines",
                    {
                        "type": "logical_group",
                        "operator": "OR",
                        "targets": [
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            "Carboxylic acid with primary amine to amide",
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                            "Acyl chloride with secondary amine to amide",
                            "Ester with primary amine to amide",
                            "Ester with secondary amine to amide",
                            "Schotten-Baumann_amide",
                            "Acylation of primary amines",
                            "Acylation of secondary amines"
                        ]
                    }
                ]
            }
        })
        # Positional constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": {
                    "type": "logical_group",
                    "operator": "OR",
                    "targets": [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines"
                    ]
                },
                "position": "first_half",
                "definition": "The reaction must occur in a step with depth <= max_depth / 2, where depth 0 is the final product."
            }
        })

    return combined_strategy, findings_json