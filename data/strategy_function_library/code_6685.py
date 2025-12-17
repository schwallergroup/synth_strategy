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

HETEROCYCLIC_RINGS = [
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "furan",
    "thiophene",
    "triazole",
    "tetrazole",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "pyrimidine",
    "pyrazine",
    "indole",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
    "purine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
]
COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Negishi coupling",
    "Stille reaction_vinyl",
    "Stille reaction_aryl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Stille reaction_vinyl OTf",
    "Stille reaction_aryl OTf",
    "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf",
    "Stille reaction_other",
    "Stille reaction_other OTf",
    "Heck terminal vinyl",
    "Heck_non-terminal_vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) coupling of two heterocyclic fragments. This is determined by checking if the reaction is a known coupling type from the `COUPLING_REACTIONS` list and if at least two of the reactants contain a ring from the `HETEROCYCLIC_RINGS` list.
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

    final_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal final_coupling_detected, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for final coupling (depth 0 or 1 - late stage)
            if depth <= 1:
                # Check if this is a coupling reaction
                is_coupling = False
                for rxn in COUPLING_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_coupling = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                
                if is_coupling and len(reactants_smiles) >= 2:
                    # Check if reactants contain heterocycles
                    heterocycle_counts = [[] for _ in range(len(reactants_smiles))]
                    reactants_with_heterocycles = 0

                    for i, reactant_smiles in enumerate(reactants_smiles):
                        # Check for predefined heterocycles
                        found_heterocycles_in_reactant = False
                        for ring in HETEROCYCLIC_RINGS:
                            if checker.check_ring(ring, reactant_smiles):
                                heterocycle_counts[i].append(ring)
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                found_heterocycles_in_reactant = True
                        if found_heterocycles_in_reactant:
                            reactants_with_heterocycles += 1

                    # If at least two reactants have heterocycles and they're being joined
                    if reactants_with_heterocycles >= 2:
                        final_coupling_detected = True
                        # Record structural constraints
                        if {"type": "positional", "details": {"target": "coupling_of_heterocyclic_fragments", "position": "late_stage (depth <= 1)"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "coupling_of_heterocyclic_fragments", "position": "late_stage (depth <= 1)"}})
                        if {"type": "count", "details": {"target": "reactants_with_heterocyclic_rings_in_coupling_step", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_with_heterocyclic_rings_in_coupling_step", "operator": ">=", "value": 2}})

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return final_coupling_detected, findings_json