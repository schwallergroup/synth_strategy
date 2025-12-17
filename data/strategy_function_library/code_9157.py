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


N_ALKYLATION_REACTIONS = [
    "N-alkylation of secondary amines with alkyl halides",
    "N-methylation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects an amide N-alkylation in the final synthetic step. The strategy is identified if a reaction is classified as a specific N-alkylation type (from the N_ALKYLATION_REACTIONS list) and it transforms a secondary amide reactant into a tertiary amide product.
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

    # Track if we found the pattern
    found_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation, findings_json

        # The final reaction step is at depth 1
        if node["type"] == "reaction" and depth == 1:
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_n_alkylation_type = False
                detected_reaction_name = None
                # Check if this is a known N-alkylation reaction type
                for rxn_name in N_ALKYLATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_n_alkylation_type = True
                        detected_reaction_name = rxn_name
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                if is_n_alkylation_type:
                    # Verify that we're alkylating an amide
                    has_sec_amide = False
                    for reactant in reactants:
                        if checker.check_fg("Secondary amide", reactant):
                            has_sec_amide = True
                            if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                            break

                    # Check product for tertiary amide
                    has_tert_amide = checker.check_fg("Tertiary amide", product)
                    if has_tert_amide:
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    if has_sec_amide and has_tert_amide:
                        found_n_alkylation = True
                        # Add structural constraints if all conditions are met
                        if detected_reaction_name == "N-alkylation of secondary amines with alkyl halides":
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "N-alkylation of secondary amines with alkyl halides",
                                        "Secondary amide",
                                        "Tertiary amide"
                                    ],
                                    "scope": "reaction_step"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "N-alkylation of secondary amines with alkyl halides",
                                    "position": "last_stage"
                                }
                            })
                        elif detected_reaction_name == "N-methylation":
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "N-methylation",
                                        "Secondary amide",
                                        "Tertiary amide"
                                    ],
                                    "scope": "reaction_step"
                                }
                            })
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "N-methylation",
                                    "position": "last_stage"
                                }
                            })

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_n_alkylation, findings_json