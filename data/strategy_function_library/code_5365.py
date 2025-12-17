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


HETEROCYCLE_COUPLING_RINGS = [
    "pyrazole", "piperazine", "morpholine", "pyridine", "pyrimidine",
    "pyrazine", "imidazole", "oxazole", "thiazole", "triazole",
    "tetrazole", "furan", "thiophene", "pyrrole", "indole",
    "quinoline", "isoquinoline", "benzimidazole", "benzoxazole",
    "benzothiazole", "pyrrolidine", "piperidine", "oxadiazole",
    "thiadiazole", "benzotriazole",
]

HETEROCYCLE_COUPLING_REACTIONS = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", "Buchwald-Hartwig",
    "Goldberg coupling", "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride", "Ullmann-Goldberg Substitution amine",
    "N-arylation_heterocycles", "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Golderg/N-arylation secondary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects convergent syntheses where two fragments, each containing a heterocycle
    from a predefined list, are coupled via a C-N bond-forming reaction from a
    specific list of named reactions (e.g., Buchwald-Hartwig, Ullmann).
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

    heterocycle_coupling_detected = False
    
    # Define the structural constraint object from the original JSON for easy access
    structural_constraint_obj = {
      "type": "co-occurrence",
      "details": {
        "targets": [
          "A reaction from the HETEROCYCLE_COUPLING_REACTIONS list",
          "At least two reactants each containing a ring from the HETEROCYCLE_COUPLING_RINGS list"
        ],
        "scope": "single_reaction_step"
      }
    }

    def dfs_traverse(node, depth):
        nonlocal heterocycle_coupling_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if we have multiple reactants (convergent)
            if len(reactants) >= 2:
                # Check for heterocycles in reactants
                heterocycle_reactants = []
                found_rings_in_reactants = set()
                for reactant in reactants:
                    for ring_name in HETEROCYCLE_COUPLING_RINGS:
                        if checker.check_ring(ring_name, reactant):
                            heterocycle_reactants.append(reactant)
                            if ring_name not in found_rings_in_reactants:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                found_rings_in_reactants.add(ring_name)
                            break

                # If we have at least 2 heterocycle-containing fragments
                if len(heterocycle_reactants) >= 2:
                    # Check if the reaction involves a specified C-N bond formation
                    for reaction_type in HETEROCYCLE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            heterocycle_coupling_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # Add the structural constraint if all conditions are met
                            if structural_constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_obj)
                            return  # Early return once detected

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical' or other non-reaction type
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal with initial depth 0
    dfs_traverse(route, 0)
    return heterocycle_coupling_detected, findings_json
