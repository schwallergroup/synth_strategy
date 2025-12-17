from typing import Tuple, Dict, List
import copy
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


# Refactoring for Enumeration: Isolate the lists of chemical entities.
HETEROCYCLES_OF_INTEREST = [
    "furan", "pyran", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine",
    "triazole", "tetrazole", "morpholine", "piperidine", "piperazine",
    "isoxazole", "isothiazole", "oxadiazole", "thiadiazole", "indole",
    "quinoline", "isoquinoline", "benzoxazole", "benzothiazole", "benzimidazole",
]

# Consolidate and clean up coupling reaction names.
COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki", "Buchwald-Hartwig", "N-arylation", "Heck", "Sonogashira",
    "Stille", "Negishi", "Kumada", "Ullmann", "Chan-Lam",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "N-arylation_heterocycles",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the late-stage coupling of two fragments where both contain a heterocycle. The reaction must be one of the specified coupling types (e.g., Suzuki, Buchwald-Hartwig) and both reactants must contain one of the specified heterocyclic rings (e.g., pyridine, indole).
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

    found_late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_coupling, findings_json

        if node.get("type") == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a late-stage reaction (depth 0 or 1)
            is_late_stage = False
            if depth <= 1:
                is_late_stage = True
                if {"type": "positional", "details": {"target": "coupling_reaction", "position": "late_stage (depth <= 1)"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "coupling_reaction", "position": "late_stage (depth <= 1)"}})

            if is_late_stage:
                # Check if we're joining two or more fragments
                is_multi_fragment = False
                if len(reactants) >= 2:
                    is_multi_fragment = True
                    if {"type": "count", "details": {"target": "reactants_in_coupling_step", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_in_coupling_step", "operator": ">=", "value": 2}})

                if is_multi_fragment:
                    # Check if this is a coupling reaction
                    is_coupling = False
                    for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_coupling = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            break

                    if is_coupling:
                        # Check if at least two reactants contain heterocycles
                        heterocycle_count = 0
                        for reactant in reactants:
                            for ring in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring, reactant):
                                    heterocycle_count += 1
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    break # Move to next reactant

                        if heterocycle_count >= 2:
                            if {"type": "count", "details": {"target": "reactants_with_heterocycle", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_with_heterocycle", "operator": ">=", "value": 2}})
                            found_late_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_late_coupling, findings_json