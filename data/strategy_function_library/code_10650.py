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


AMINE_FUNCTIONALIZATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "Reductive amination with alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves late-stage amine functionalization, specifically the conversion of a primary or secondary amine to a tertiary amine. This is identified by checking for one of the following reaction types in the final steps: N-alkylation of primary amines with alkyl halides, N-alkylation of secondary amines with alkyl halides, Methylation with MeI_primary, Methylation with MeI_secondary, Reductive amination with aldehyde, Reductive amination with ketone, Eschweiler-Clarke Primary Amine Methylation, Eschweiler-Clarke Secondary Amine Methylation, Reductive methylation of primary amine with formaldehyde, and Reductive amination with alcohol.
    """
    amine_functionalization_found = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal amine_functionalization_found, findings_json
        if amine_functionalization_found:
            return

        # Consider depths 0-2 as late-stage
        is_late_stage = depth <= 2

        if (
            node["type"] == "reaction"
            and is_late_stage
            and "metadata" in node
            and "mapped_reaction_smiles" in node["metadata"]
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for specific amine functionalization reactions
            is_amine_functionalization = False
            for reaction_name in AMINE_FUNCTIONALIZATION_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_amine_functionalization = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)

            if is_amine_functionalization:
                # Confirm conversion from primary/secondary to tertiary amine
                reactant_has_pri_or_sec_amine = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r):
                        reactant_has_pri_or_sec_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", r):
                        reactant_has_pri_or_sec_amine = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                product_has_tertiary_amine = checker.check_fg("Tertiary amine", product)
                if product_has_tertiary_amine:
                    if "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")

                if reactant_has_pri_or_sec_amine and product_has_tertiary_amine:
                    amine_functionalization_found = True
                    # Add structural constraint for co-occurrence
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "any_amine_functionalization_reaction",
                                "reactant_with_primary_or_secondary_amine",
                                "product_with_tertiary_amine"
                            ],
                            "scope": "reaction_step"
                        }
                    })
                    # Add structural constraint for positional (late_stage)
                    if is_late_stage:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "amine_functionalization_and_conversion_event",
                                "position": "late_stage"
                            }
                        })
                    return

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amine_functionalization_found, findings_json
