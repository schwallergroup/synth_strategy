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


HETEROCYCLES_OF_INTEREST = [
    "pyrazole", "imidazole", "oxazole", "thiazole", "triazole", "tetrazole",
    "isoxazole", "isothiazole", "oxadiazole", "thiadiazole", "benzoxazole",
    "benzothiazole", "benzimidazole", "indole", "pyrrole", "furan", "thiophene",
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "quinoline", "isoquinoline",
    "piperidine", "piperazine", "morpholine", "thiomorpholine",
]

HETEROCYCLE_REACTION_MAP = {
    "pyrazole": ["pyrazole", "{pyrazole}", "Pyrazole formation"],
    "imidazole": ["imidazole", "{imidazole}", "{triaryl-imidazole}"],
    "oxazole": ["oxazole"],
    "thiazole": ["thiazole", "{thiazole}"],
    "triazole": [
        "triazole", "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition", "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{Huisgen_Ru-catalyzed_1,5_subst}",
    ],
    "tetrazole": [
        "tetrazole", "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "Azide-nitrile click cycloaddition to tetrazole",
    ],
    "isoxazole": ["isoxazole"],
    "isothiazole": ["isothiazole"],
    "oxadiazole": [
        "oxadiazole", "{oxadiazole}",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    ],
    "thiadiazole": ["thiadiazole"],
    "benzoxazole": [
        "benzoxazole", "{benzoxazole}", "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)", "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
    ],
    "benzothiazole": [
        "benzothiazole", "{benzothiazole}", "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
    ],
    "benzimidazole": [
        "benzimidazole", "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}", "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid", "Benzimidazole aldehyde",
    ],
    "indole": [
        "indole", "{indole}", "{Fischer indole}", "Paal-Knorr pyrrole synthesis",
        "{piperidine_indole}",
    ],
    "pyrrole": ["pyrrole", "{Paal-Knorr pyrrole}", "Paal-Knorr pyrrole synthesis"],
    "furan": ["furan", "{benzofuran}"],
    "thiophene": ["thiophene", "{benzothiophene}"],
}

# Add missing heterocycles to the reactions dictionary
for hc in HETEROCYCLES_OF_INTEREST:
    if hc not in HETEROCYCLE_REACTION_MAP:
        HETEROCYCLE_REACTION_MAP[hc] = [hc]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if a heterocycle is formed
    in the late stage of the synthesis.
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

    # Track formation depths for each heterocycle
    heterocycle_formation_depths = {hc: -1 for hc in HETEROCYCLES_OF_INTEREST}
    max_depth = -1
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycle formation
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    # Check if heterocycle exists in reactants and product
                    has_heterocycle_in_reactants = any(
                        checker.check_ring(heterocycle, r) for r in reactants
                    )
                    if has_heterocycle_in_reactants:
                        # Add to findings if found in reactants (even if not formed)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    has_heterocycle_in_product = checker.check_ring(heterocycle, product)
                    if has_heterocycle_in_product:
                        # Add to findings if found in product
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    # Detect heterocycle formation (not in reactants but in product)
                    if not has_heterocycle_in_reactants and has_heterocycle_in_product:
                        # Check for specific heterocycle formation reactions
                        if heterocycle in HETEROCYCLE_REACTION_MAP:
                            for reaction_name in HETEROCYCLE_REACTION_MAP[heterocycle]:
                                if checker.check_reaction(reaction_name, rsmi):
                                    heterocycle_formation_depths[heterocycle] = depth
                                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                                    break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if any heterocycle formation is in the late stage of the synthesis
    late_stage_threshold = max_depth / 2

    for heterocycle, depth in heterocycle_formation_depths.items():
        if depth != -1:
            is_late_stage = depth <= late_stage_threshold
            if is_late_stage:
                result = True
                # Add the structural constraint if a late-stage heterocycle formation is found
                if {"type": "positional", "details": {"target": "heterocycle_formation", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "heterocycle_formation",
                            "position": "late_stage"
                        }
                    })

    return result, findings_json
