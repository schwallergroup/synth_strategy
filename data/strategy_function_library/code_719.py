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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves multiple functional groups common in thiol manipulation pathways
    (e.g., alcohol, mesylate, thiol, sulfide, sulfone) and contains at least one characteristic reaction
    like S-alkylation or sulfur oxidation.
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

    # Track the functional group states we observe
    observed_states = {
        "alcohol": False,  # Primary, Secondary, or Tertiary alcohol
        "mesylate": False,  # Mesylate
        "thioacetate": False,  # Thioacetate
        "thiol": False,  # Aliphatic or Aromatic thiol
        "methylthioether": False,  # Monosulfide
        "sulfone": False,  # Sulfone
    }

    # Track if we've seen relevant thiol manipulation reactions
    thiol_reactions_observed = False
    observed_thiol_reactions_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal thiol_reactions_observed, observed_thiol_reactions_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check reactants for functional groups
            for r in reactants:
                if checker.check_fg("Primary alcohol", r):
                    observed_states["alcohol"] = True
                    if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                if checker.check_fg("Secondary alcohol", r):
                    observed_states["alcohol"] = True
                    if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                if checker.check_fg("Tertiary alcohol", r):
                    observed_states["alcohol"] = True
                    if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                if checker.check_fg("Mesylate", r):
                    observed_states["mesylate"] = True
                    if "Mesylate" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
                if checker.check_fg("Aliphatic thiol", r):
                    observed_states["thiol"] = True
                    if "Aliphatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aliphatic thiol")
                if checker.check_fg("Aromatic thiol", r):
                    observed_states["thiol"] = True
                    if "Aromatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aromatic thiol")
                if checker.check_fg("Monosulfide", r):
                    observed_states["methylthioether"] = True
                    if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
                if checker.check_fg("Sulfone", r):
                    observed_states["sulfone"] = True
                    if "Sulfone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfone")

            # Check product for functional groups
            if checker.check_fg("Primary alcohol", product):
                observed_states["alcohol"] = True
                if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
            if checker.check_fg("Secondary alcohol", product):
                observed_states["alcohol"] = True
                if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
            if checker.check_fg("Tertiary alcohol", product):
                observed_states["alcohol"] = True
                if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
            if checker.check_fg("Mesylate", product):
                observed_states["mesylate"] = True
                if "Mesylate" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
            if checker.check_fg("Aliphatic thiol", product):
                observed_states["thiol"] = True
                if "Aliphatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aliphatic thiol")
            if checker.check_fg("Aromatic thiol", product):
                observed_states["thiol"] = True
                if "Aromatic thiol" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic thiol")
            if checker.check_fg("Monosulfide", product):
                observed_states["methylthioether"] = True
                if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")
            if checker.check_fg("Sulfone", product):
                observed_states["sulfone"] = True
                if "Sulfone" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfone")

            # Check for specific thiol-related reactions
            thiol_reaction_names = [
                "S-alkylation of thiols",
                "S-alkylation of thiols (ethyl)",
                "S-alkylation of thiols with alcohols",
                "S-alkylation of thiols with alcohols (ethyl)",
                "Sulfanyl to sulfinyl",
                "Sulfanyl to sulfinyl_peroxide",
                "Sulfanyl to sulfinyl_H2O2",
                "Sulfanyl to sulfinyl_SO3-",
                "Sulfanyl to sulfinyl_sulfonyl",
                "Formation of Sulfonic Esters"
            ]
            for r_name in thiol_reaction_names:
                if checker.check_reaction(r_name, rsmi):
                    thiol_reactions_observed = True
                    observed_thiol_reactions_count += 1
                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we've observed at least 3 of the states AND relevant reactions
    states_count = sum(1 for state in observed_states.values() if state)

    result = states_count >= 3 and thiol_reactions_observed

    # Record structural constraints if met
    if states_count >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_thiol_pathway_functional_group_classes",
                "operator": ">=",
                "value": 3
            }
        })
    if thiol_reactions_observed:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "thiol_manipulation_reactions",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json
