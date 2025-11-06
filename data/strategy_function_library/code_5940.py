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

HETEROCYCLES_OF_INTEREST = [
    "pyrazole", "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "furan", "thiophene", "oxazole", "thiazole", "isoxazole",
    "isothiazole", "imidazole", "oxadiazole", "thiadiazole", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole", "quinoline", "isoquinoline",
    "purine", "piperazine", "morpholine", "piperidine", "pyrrolidine"
]

FRAGMENT_COUPLING_REACTIONS = [
    "Suzuki", "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf", "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters", "Negishi", "Stille", "Stille reaction_vinyl",
    "Stille reaction_aryl", "Stille reaction_benzyl", "Stille reaction_allyl",
    "Stille reaction_vinyl OTf", "Stille reaction_aryl OTf", "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf", "Stille reaction_other", "Stille reaction_other OTf",
    "Heck", "Heck terminal vinyl", "Heck_terminal_vinyl", "Heck_non-terminal_vinyl",
    "Oxidative Heck reaction", "Oxidative Heck reaction with vinyl ester",
    "Heck reaction with vinyl ester and amine", "Sonogashira",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf", "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_alkenyl halide", "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl OTf", "Sonogashira alkyne_alkenyl OTf",
    "Sonogashira acetylene_acyl halide", "Sonogashira alkyne_acyl halide",
    "Buchwald-Hartwig", "N-arylation", "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", "Ullmann-Goldberg",
    "Ullmann-Goldberg Substitution amine", "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol", "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride", "Goldberg coupling aryl amide-aryl chloride",
    "Ullmann condensation", "decarboxylative_coupling", "Hiyama-Denmark Coupling",
    "Kumada cross-coupling", "Aryllithium cross-coupling", "Ackermann Reaction",
    "Catellani reaction ortho", "Catellani reaction para", "beta C(sp3) arylation",
    "C-methylation", "Chan-Lam alcohol", "Chan-Lam amine", "Chan-Lam etherification"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) fragment coupling reaction where at least two reactant fragments each contain a heterocycle. The strategy is confirmed if the reaction is a known cross-coupling type from a predefined list (e.g., Suzuki, Stille, Buchwald-Hartwig) and the reactants contain heterocycles from a specific list (e.g., pyridine, pyrazole, indole).
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

    heterocycle_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_coupling, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "heterocycle_fragment_coupling",
                    "position": "late_stage",
                    "max_depth": 2
                }
            })
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) >= 2:  # Need at least two reactants for coupling
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants_in_coupling_reaction",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    heterocycle_reactants_count = 0
                    for reactant in reactants:
                        reactant_has_heterocycle = False
                        for ring_name in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring_name, reactant):
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                heterocycle_reactants_count += 1
                                reactant_has_heterocycle = True
                                break

                    if heterocycle_reactants_count >= 2:
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "reactants_with_target_heterocycle",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        # Check if this is actually a coupling reaction
                        for rxn_type in FRAGMENT_COUPLING_REACTIONS:
                            if checker.check_reaction(rxn_type, rsmi):
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                heterocycle_coupling = True
                                return

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return heterocycle_coupling, findings_json