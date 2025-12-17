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


NARYLATION_REACTION_TYPES = [
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "N-arylation_heterocycles",
    "Buchwald-Hartwig",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
]

HETEROCYCLE_FORMATION_REACTION_TYPES = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "imidazole",
    "pyrazole",
    "benzimidazole formation from aldehyde",
    "benzimidazole formation from acyl halide",
    "benzimidazole formation from ester/carboxylic acid",
    "benzoxazole formation from aldehyde",
    "benzoxazole formation from acyl halide",
    "benzoxazole formation from ester/carboxylic acid",
    "benzoxazole formation (intramolecular)",
    "benzothiazole formation from aldehyde",
    "benzothiazole formation from acyl halide",
    "benzothiazole formation from ester/carboxylic acid",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a heterocycle is formed in an early stage, followed by a late-stage N-arylation for fragment coupling. This function uses the curated lists `HETEROCYCLE_FORMATION_REACTION_TYPES` and `NARYLATION_REACTION_TYPES` to identify the respective reaction steps and compares their positions in the synthetic sequence.
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

    heterocycle_formation_depth = -1
    narylation_depths = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for heterocycle formation reactions
            for rxn in HETEROCYCLE_FORMATION_REACTION_TYPES:
                if checker.check_reaction(rxn, rsmi):
                    if heterocycle_formation_depth < depth:
                        heterocycle_formation_depth = depth
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)

            # Check for N-arylation using checker functions
            for rxn in NARYLATION_REACTION_TYPES:
                if checker.check_reaction(rxn, rsmi):
                    narylation_depths.append(depth)
                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)

        # Traverse children (deeper in retrosynthesis = earlier in forward synthesis)
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Strategy requires heterocycle formation to be early-stage (higher depth)
    # and N-arylation to be late-stage (lower depth)
    if heterocycle_formation_depth > 0:
        # Check if at least one N-arylation occurs at a lower depth (later stage)
        # than the heterocycle formation (early stage)
        for narylation_depth in narylation_depths:
            if narylation_depth < heterocycle_formation_depth:
                result = True
                # Add the structural constraint if the condition is met
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before_any_of": [
                            "Formation of NOS Heterocycles",
                            "Paal-Knorr pyrrole synthesis",
                            "benzimidazole_derivatives_carboxylic-acid/ester",
                            "benzimidazole_derivatives_aldehyde",
                            "benzothiazole",
                            "benzoxazole_arom-aldehyde",
                            "benzoxazole_carboxylic-acid",
                            "thiazole",
                            "tetrazole_terminal",
                            "imidazole",
                            "pyrazole",
                            "benzimidazole formation from aldehyde",
                            "benzimidazole formation from acyl halide",
                            "benzimidazole formation from ester/carboxylic acid",
                            "benzoxazole formation from aldehyde",
                            "benzoxazole formation from acyl halide",
                            "benzoxazole formation from ester/carboxylic acid",
                            "benzoxazole formation (intramolecular)",
                            "benzothiazole formation from aldehyde",
                            "benzothiazole formation from acyl halide",
                            "benzothiazole formation from ester/carboxylic acid"
                        ],
                        "after_any_of": [
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "N-arylation_heterocycles",
                            "Buchwald-Hartwig",
                            "Ullmann-Goldberg Substitution amine",
                            "Goldberg coupling",
                            "Goldberg coupling aryl amine-aryl chloride",
                            "Goldberg coupling aryl amide-aryl chloride"
                        ]
                    }
                })
                break # Only need one instance to satisfy the condition

    return result, findings_json
