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

# Module-level constants for better organization and clarity
HETEROCYCLE_TYPES = [
    "triazole", "tetrazole", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "indole", "benzimidazole",
    "benzoxazole", "benzothiazole", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "benzotriazole",
]

FUNCTIONALIZATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of primary amines", "Acylation of secondary amines",
    "Methylation with MeI_primary", "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary", "N-methylation", "O-methylation",
    "Schotten-Baumann to ester",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Esterification of Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
]

COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Sonogashira alkyne_aryl halide", "Heck terminal vinyl",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Stille reaction_aryl", "Negishi coupling", "Hiyama-Denmark Coupling",
    "Kumada cross-coupling", "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a 'build-couple-pair' strategy where a complex heterocycle is built
    before final functionalization steps.

    In retrosynthetic analysis:
    - Higher depth = earlier stage in forward synthesis (build phase)
    - Lower depth = later stage in forward synthesis (pair/functionalization phase)
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

    # We need to track multiple aspects of the synthesis
    heterocycle_built = False
    coupling_phase = False
    late_functionalization = False
    heterocycle_build_depth = -1
    coupling_depth = -1
    functionalization_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_built, coupling_phase, late_functionalization
        nonlocal heterocycle_build_depth, coupling_depth, functionalization_depth
        nonlocal findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycle formation (build phase)
                if not heterocycle_built:
                    product_has_heterocycle = False
                    for heterocycle in HETEROCYCLE_TYPES:
                        if checker.check_ring(heterocycle, product):
                            product_has_heterocycle = True
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break

                    if product_has_heterocycle:
                        reactants_have_heterocycle = False
                        for reactant in reactants:
                            for heterocycle in HETEROCYCLE_TYPES:
                                if checker.check_ring(heterocycle, reactant):
                                    reactants_have_heterocycle = True
                                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                    break
                            if reactants_have_heterocycle:
                                break

                        if not reactants_have_heterocycle:
                            heterocycle_built = True
                            heterocycle_build_depth = depth
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check for coupling phase (couple phase)
                if not coupling_phase:
                    for reaction in COUPLING_REACTIONS:
                        if checker.check_reaction(reaction, rsmi):
                            coupling_phase = True
                            coupling_depth = depth
                            if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction)
                            break

                # Check for late functionalization (pair phase)
                if not late_functionalization:
                    for reaction in FUNCTIONALIZATION_REACTIONS:
                        if checker.check_reaction(reaction, rsmi):
                            late_functionalization = True
                            functionalization_depth = depth
                            if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction)
                            break

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # For a true build-couple-pair strategy:
    # 1. Heterocycle building must occur before functionalization.
    # 2. In retrosynthetic analysis, this means build_depth > functionalization_depth.
    valid_strategy = False
    if heterocycle_built and late_functionalization:
        # Add the co-occurrence constraint if both are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ring_formation",
                    "functionalization_reaction"
                ],
                "description": "The strategy requires both the formation of a heterocycle (from HETEROCYCLE_TYPES) and a late-stage functionalization reaction (from FUNCTIONALIZATION_REACTIONS) to occur."
            }
        })
        if heterocycle_build_depth > functionalization_depth:
            valid_strategy = True

    # Coupling is optional but should be between build and pair phases if present
    if coupling_phase and valid_strategy:
        if not (functionalization_depth < coupling_depth < heterocycle_build_depth):
            # If coupling is present but out of order, the overall strategy is invalid
            valid_strategy = False

    # Add the sequence constraint if the strategy is valid
    if valid_strategy:
        sequence_constraint = {
            "type": "sequence",
            "details": {
                "description": "The heterocycle formation must occur at a greater depth (earlier in synthesis) than the functionalization. If a coupling reaction occurs, it must happen between the ring formation and the functionalization. The depth check is `functionalization_depth < coupling_depth < heterocycle_build_depth`.",
                "ordered_events": [
                    "ring_formation",
                    "coupling_reaction",
                    "functionalization_reaction"
                ],
                "optional_events": [
                    "coupling_reaction"
                ]
            }
        }
        # Only add if not already present (e.g., from a previous check if logic was different)
        if sequence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(sequence_constraint)

    return valid_strategy, findings_json
