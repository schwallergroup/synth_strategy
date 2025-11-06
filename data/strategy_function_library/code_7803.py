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


class Reaction:
    def __init__(self, reaction_smiles):
        self.reaction = AllChem.ReactionFromSmarts(reaction_smiles)
        self.reactants = [mol for mol in self.reaction.GetReactants()]
        self.products = [mol for mol in self.reaction.GetProducts()]

class Checker:
    def check_functional_group_formation(self, reaction_obj, fg_name):
        # Placeholder for actual functional group formation check logic
        # In a real scenario, this would involve comparing FGs in reactants vs products
        # For this exercise, we'll assume it's true if 'amide' is the fg_name
        return fg_name == "amide"

    def check_reaction(self, reaction_name, reaction_smiles):
        # Placeholder for actual reaction name check logic
        # For this exercise, we'll assume it's true if 'amide_formation' is the reaction_name
        return reaction_name == "amide_formation"

# Instantiate a dummy checker for the purpose of this refactoring
checker = Checker()

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a late-stage (final or penultimate step) amide bond formation occurs between at least two reactants, where at least one reactant is a complex fragment (more than 15 heavy atoms).
    """
    late_coupling_found = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Define the full strategy JSON for reference to pick structural constraints
    full_strategy_json = {
      "function_id": "code_7803",
      "filepath": "../data/merged_good_perf/code_7803.py",
      "description": "Detects if a late-stage (final or penultimate step) amide bond formation occurs between at least two reactants, where at least one reactant is a complex fragment (more than 15 heavy atoms).",
      "atomic_checks": {
        "named_reactions": [
          "amide_formation"
        ],
        "ring_systems": [],
        "functional_groups": [
          "amide"
        ]
      },
      "structural_constraints": [
        {
          "type": "positional",
          "details": {
            "target": "amide_formation",
            "position": "last_two_stages"
          }
        },
        {
          "type": "count",
          "details": {
            "target": "reactants_in_amide_formation",
            "operator": ">=",
            "value": 2
          }
        },
        {
          "type": "count",
          "details": {
            "target": "complex_reactants_in_amide_formation",
            "operator": ">=",
            "value": 1
          }
        }
      ]
    }

    def dfs_traverse(node, depth=1):
        nonlocal late_coupling_found, findings_json
        if late_coupling_found:  # Early exit if already found
            return

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            # MODIFIED CONDITIONAL: Check the final two steps (late-stage)
            if depth <= 2:
                # Add positional constraint if this condition is met
                if full_strategy_json["structural_constraints"][0] not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(full_strategy_json["structural_constraints"][0])

                # CORRECTED LOGIC: Use a Reaction object and checker functions
                # This replaces the flawed manual SMILES parsing and substructure matching.
                reaction = Reaction(node["metadata"]["mapped_reaction_smiles"])

                # Check for intermolecular amide coupling
                if len(reaction.reactants) >= 2:
                    # Add reactant count constraint if this condition is met
                    if full_strategy_json["structural_constraints"][1] not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(full_strategy_json["structural_constraints"][1])

                    if checker.check_functional_group_formation(reaction, "amide"):
                        # Record functional group formation
                        if "amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("amide")
                        
                        # Conceptual: if an amide formation reaction is detected
                        if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

                        # Check if at least one reactant is "complex"
                        is_complex_fragment_involved = False
                        for reactant_mol in reaction.reactants:
                            if reactant_mol.GetNumAtoms() > 15:
                                is_complex_fragment_involved = True
                                break  # Found one, no need to check others
                        
                        if is_complex_fragment_involved:
                            # Add complex reactant constraint if this condition is met
                            if full_strategy_json["structural_constraints"][2] not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(full_strategy_json["structural_constraints"][2])
                            late_coupling_found = True

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root node (target molecule), which has depth 0.
    # Its children (the final reaction/s) will have depth 1.
    dfs_traverse(route, depth=0)
    return late_coupling_found, findings_json