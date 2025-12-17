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


# Refactoring for Enumeration: Define lists at the module level
NITROGEN_HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrrole", "imidazole", "pyrazole", "triazole", "tetrazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "indole",
    "quinoline", "isoquinoline", "benzimidazole", "benzotriazole",
]

CN_COUPLING_REACTIONS_OF_INTEREST = [
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "N-arylation_heterocycles",
    "Goldberg coupling",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage C-N coupling between an aryl halide and a nitrogen heterocycle.
    The reaction must match a name from the `CN_COUPLING_REACTIONS_OF_INTEREST` list,
    and the heterocycle must be present in the `NITROGEN_HETEROCYCLES_OF_INTEREST` list.
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

    c_n_coupling_detected = False

    # Define late stage depth threshold
    LATE_STAGE_DEPTH = 2  # Consider depths 1 and 2 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal c_n_coupling_detected, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction
            if depth <= LATE_STAGE_DEPTH:
                # Add positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "C-N coupling reaction",
                        "position_type": "depth",
                        "operator": "<=",
                        "value": 2
                    }
                })

                # Find reactants with aryl halides
                aryl_halide_reactants = []
                for r in reactants:
                    if r and checker.check_fg("Aromatic halide", r):
                        aryl_halide_reactants.append(r)
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                # Find reactants with nitrogen heterocycles
                n_heterocycle_reactants = []
                found_n_heterocycle_ring = False
                for r in reactants:
                    for ring in NITROGEN_HETEROCYCLES_OF_INTEREST:
                        if r and checker.check_ring(ring, r):
                            n_heterocycle_reactants.append(r)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            found_n_heterocycle_ring = True
                            break # Found one ring for this reactant, move to next reactant

                # Check if we have both required components
                has_aryl_halide = len(aryl_halide_reactants) > 0
                has_n_heterocycle = len(n_heterocycle_reactants) > 0

                # First check if this is a known C-N coupling reaction
                is_cn_coupling = False
                for rxn in CN_COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn, rsmi):
                        is_cn_coupling = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                # If we have all components and it's a C-N coupling reaction
                if has_aryl_halide and has_n_heterocycle and is_cn_coupling:
                    # Add co-occurrence constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "scope": "reaction_step",
                            "targets": [
                                "named_cn_coupling",
                                "reactant_with_aromatic_halide",
                                "reactant_with_nitrogen_heterocycle"
                            ]
                        }
                    })

                    # A successful C-N coupling must consume the aryl halide.
                    if not checker.check_fg("Aromatic halide", product):
                        c_n_coupling_detected = True
                        # Add negation constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "Aromatic halide",
                                "scope": "product"
                            }
                        })

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return c_n_coupling_detected, findings_json
