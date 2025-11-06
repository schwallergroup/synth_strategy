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


# Refactored lists for enumeration
EARLY_STAGE_RING_FORMATION_REACTIONS = [
    "Paal-Knorr pyrrole synthesis", "Fischer indole", "Friedlaender chinoline",
    "benzofuran", "benzothiophene", "indole", "Diels-Alder", "Diels-Alder (ON bond)",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
    "benzimidazole_derivatives_aldehyde", "benzothiazole", "benzoxazole_arom-aldehyde",
    "thiazole", "Pictet-Spengler", "Niementowski_quinazoline", "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide", "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine", "spiro-chromanone", "pyrazole", "phthalazinone",
    "triaryl-imidazole", "oxadiazole", "imidazole", "Pauson-Khand reaction",
    "Azide-nitrile click cycloaddition to tetrazole", "Azide-nitrile click cycloaddition to triazole",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    "Intramolecular transesterification/Lactone formation", "Formation of NOS Heterocycles",
]

MID_STAGE_BIARYL_COUPLING_REACTIONS = [
    "Suzuki", "Suzuki coupling with boronic acids", "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters", "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with boronic esters", "Stille", "Stille reaction_aryl",
    "Stille reaction_aryl OTf", "Negishi", "Kumada cross-coupling", "Hiyama-Denmark Coupling",
    "Ullmann condensation", "Ullmann-Goldberg Substitution aryl alcohol",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", "Buchwald-Hartwig",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride", "Goldberg coupling aryl amide-aryl chloride",
    "Aryllithium cross-coupling",
]

LATE_STAGE_FG_MODIFICATION_REACTIONS = [
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Oxidation of aldehydes to carboxylic acids",
    "Reduction of aldehydes and ketones to alcohols", "Reductive amination with aldehyde",
    "Reductive amination with ketone", "Acylation of primary amines",
    "Acylation of secondary amines",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of alcohol to carboxylic acid", "Oxidation of ketone to carboxylic acid",
    "Reduction of ester to primary alcohol", "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol", "Nitrile to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide", "Ester with primary amine to amide",
    "Ester saponification (methyl deprotection)", "Ester saponification (alkyl deprotection)",
    "Methylation", "N-methylation", "O-methylation", "Alcohol to chloride",
    "Alcohol to triflate conversion", "Carboxylic acid to amide conversion",
    "Schotten-Baumann to ester", "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Boc amine deprotection", "Boc amine protection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis involves an early-stage ring formation (at high depth)
    followed by a mid-synthesis biaryl coupling and late-stage functional group modifications.
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

    ring_formation_depth = -1
    biaryl_coupling_depth = -1
    late_fg_mod = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)

    # A multi-stage strategy is not meaningful for very short routes
    if max_depth < 4:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "route_depth",
                "operator": ">=",
                "value": 4
            }
        })
        return False, findings_json
    else:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "route_depth",
                "operator": ">=",
                "value": 4
            }
        })

    # Define early, mid, and late stage based on max_depth
    early_stage_threshold = max(max_depth - 3, 2)
    mid_stage_threshold = max(max_depth // 2, 1)
    late_stage_threshold = 1

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_depth, biaryl_coupling_depth, late_fg_mod, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for ring formation at early stage (only find the first one)
            if depth >= early_stage_threshold and ring_formation_depth == -1:
                for rxn_name in EARLY_STAGE_RING_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        ring_formation_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

            # Check for biaryl coupling at mid-stage (only find the first one)
            if mid_stage_threshold - 1 <= depth <= mid_stage_threshold + 1 and biaryl_coupling_depth == -1:
                for rxn_name in MID_STAGE_BIARYL_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        biaryl_coupling_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

            # Check for late-stage functional group modifications (flag if any)
            if depth <= late_stage_threshold and not late_fg_mod:
                for rxn_name in LATE_STAGE_FG_MODIFICATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        late_fg_mod = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if all components of the strategy were detected in their respective stages
    found_early_ring_formation = ring_formation_depth != -1
    found_mid_biaryl_coupling = biaryl_coupling_depth != -1
    
    result = found_early_ring_formation and found_mid_biaryl_coupling and late_fg_mod

    if found_early_ring_formation:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    {
                        "name": "early_stage_ring_formation",
                        "position_constraint": {
                            "position": "early_stage",
                            "definition": "depth >= max(max_depth - 3, 2)"
                        },
                        "reaction_group": []
                    }
                ]
            }
        })
        # Add the specific reaction name if found, to the reaction_group of the constraint
        # This part is tricky as the original JSON doesn't specify which reaction was found.
        # For now, we'll just add the constraint itself.

    if found_mid_biaryl_coupling:
        # If the sequence constraint already exists, append to it, otherwise create a new one.
        # This assumes the sequence constraint is built incrementally.
        # For simplicity, we'll add a new one for each found component, which might duplicate the 'type: sequence' entry
        # but correctly reflects the found components.
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    {
                        "name": "mid_stage_biaryl_coupling",
                        "position_constraint": {
                            "position": "mid_stage",
                            "definition": "depth is between (max(max_depth / 2, 1) - 1) and (max(max_depth / 2, 1) + 1)"
                        },
                        "reaction_group": []
                    }
                ]
            }
        })

    if late_fg_mod:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "ordered_events": [
                    {
                        "name": "late_stage_fg_modification",
                        "position_constraint": {
                            "position": "late_stage",
                            "definition": "depth <= 1"
                        },
                        "reaction_group": []
                    }
                ]
            }
        })

    # To avoid duplicate 'type: sequence' entries, we can consolidate them.
    # This is a more robust way to handle the structural constraints.
    consolidated_sequence_constraint = {
        "type": "sequence",
        "details": {
            "ordered_events": []
        }
    }

    if found_early_ring_formation:
        consolidated_sequence_constraint["details"]["ordered_events"].append({
            "name": "early_stage_ring_formation",
            "position_constraint": {
                "position": "early_stage",
                "definition": "depth >= max(max_depth - 3, 2)"
            },
            "reaction_group": []
        })
    if found_mid_biaryl_coupling:
        consolidated_sequence_constraint["details"]["ordered_events"].append({
            "name": "mid_stage_biaryl_coupling",
            "position_constraint": {
                "position": "mid_stage",
                "definition": "depth is between (max(max_depth / 2, 1) - 1) and (max(max_depth / 2, 1) + 1)"
            },
            "reaction_group": []
        })
    if late_fg_mod:
        consolidated_sequence_constraint["details"]["ordered_events"].append({
            "name": "late_stage_fg_modification",
            "position_constraint": {
                "position": "late_stage",
                "definition": "depth <= 1"
            },
            "reaction_group": []
        })

    # Remove any existing sequence constraints and add the consolidated one if events were found
    findings_json["structural_constraints"] = [c for c in findings_json["structural_constraints"] if c.get("type") != "sequence"]
    if consolidated_sequence_constraint["details"]["ordered_events"]:
        findings_json["structural_constraints"].append(consolidated_sequence_constraint)

    return result, findings_json
