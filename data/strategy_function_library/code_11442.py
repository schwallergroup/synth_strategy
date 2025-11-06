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

# Module-level constants for enumeration
HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyrazole", "triazole", "tetrazole", "pyridine", "pyrimidine", "pyrazine",
    "pyridazine", "imidazole", "oxazole", "thiazole", "furan", "thiophene",
    "pyrrole", "indole", "benzimidazole",
]

CROSS_COUPLING_RXN_TYPES = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters", "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with boronic esters", "Negishi coupling", "Stille reaction_vinyl",
    "Stille reaction_aryl", "Stille reaction_benzyl", "Stille reaction_allyl",
    "Hiyama-Denmark Coupling", "Buchwald-Hartwig", "Ullmann condensation",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "{Suzuki}", "{N-arylation_heterocycles}", "{Buchwald-Hartwig}",
]

LATE_STAGE_RXN_TYPES = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Goldberg coupling aryl amine-aryl chloride", "Goldberg coupling aryl amide-aryl chloride",
    "Goldberg coupling", "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol", "Ullmann-Goldberg Substitution aryl alcohol",
    "Amidation reactions", "Acylation of primary amines", "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy involving the construction of specific heterocyclic systems through at least two sequential cross-coupling reactions. The strategy requires the presence of boronate intermediates (acids or esters) at some point in the route. The specific heterocycles and cross-coupling reactions checked are defined in the HETEROCYCLIC_RINGS_OF_INTEREST and CROSS_COUPLING_RXN_TYPES lists, respectively.
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

    # Track cross-coupling reactions and their depths
    cross_coupling_reactions = []
    boronate_intermediates_found = False
    late_stage_functionalization = False
    heterocyclic_systems = False

    def dfs_traverse(node, depth=0):
        nonlocal cross_coupling_reactions, boronate_intermediates_found
        nonlocal late_stage_functionalization, heterocyclic_systems, findings_json

        if node["type"] == "mol":
            # Check for heterocyclic systems in molecules
            mol_smiles = node["smiles"]
            for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                if checker.check_ring(ring, mol_smiles):
                    heterocyclic_systems = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)

            # Check for boronate intermediates
            if checker.check_fg("Boronic acid", mol_smiles):
                boronate_intermediates_found = True
                if "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
            if checker.check_fg("Boronic ester", mol_smiles):
                boronate_intermediates_found = True
                if "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for cross-coupling reactions
                for rxn_type in CROSS_COUPLING_RXN_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        cross_coupling_reactions.append((depth, rxn_type))
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                # Check for late-stage functionalization (depth 0, 1, or 2)
                if depth <= 2:
                    for rxn_type in LATE_STAGE_RXN_TYPES:
                        if checker.check_reaction(rxn_type, rsmi):
                            late_stage_functionalization = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            break

        # Traverse children with modified depth calculation
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if cross-coupling reactions are sequential
    sequential_cross_coupling = False
    if len(cross_coupling_reactions) >= 2:
        # Add structural constraint for count of cross-coupling reactions
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "cross_coupling_reaction",
                "operator": ">=",
                "value": 2,
                "description": "Requires at least two reactions from the CROSS_COUPLING_RXN_TYPES list."
            }
        })

        cross_coupling_reactions.sort(key=lambda x: x[0])
        for i in range(len(cross_coupling_reactions) - 1):
            depth_diff = abs(cross_coupling_reactions[i + 1][0] - cross_coupling_reactions[i][0])
            if depth_diff <= 3:
                sequential_cross_coupling = True
                # Add structural constraint for sequential cross-coupling
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "targets": [
                            "cross_coupling_reaction",
                            "cross_coupling_reaction"
                        ],
                        "max_depth_difference": 3,
                        "description": "At least two cross-coupling reactions must occur within 3 steps of each other."
                    }
                })
                break

    # Determine if the strategy is present
    strategy_present = (
        len(cross_coupling_reactions) >= 2
        and sequential_cross_coupling
        and boronate_intermediates_found
        and heterocyclic_systems
    )

    # Add structural constraint for boronate intermediates if found
    if boronate_intermediates_found:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain at least one molecule with a 'Boronic acid' or 'Boronic ester' functional group.",
                "any_of": [
                    "Boronic acid",
                    "Boronic ester"
                ]
            }
        })

    # Add structural constraint for heterocyclic systems if found
    if heterocyclic_systems:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain at least one molecule with a ring from the HETEROCYCLIC_RINGS_OF_INTEREST list.",
                "any_of": [
                    "pyrazole",
                    "triazole",
                    "tetrazole",
                    "pyridine",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "furan",
                    "thiophene",
                    "pyrrole",
                    "indole",
                    "benzimidazole"
                ]
            }
        })

    return strategy_present, findings_json
