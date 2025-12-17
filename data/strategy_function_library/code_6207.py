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


AROMATIC_CN_COUPLING_REACTIONS = [
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Goldberg coupling",
    "Ullmann-Goldberg Substitution amine",
    "Goldberg coupling aryl amine-aryl chloride",
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
]

HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "thiophene", "furan", "pyrrole", "oxazole", "thiazole",
    "imidazole", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
    "quinoline", "isoquinoline", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage aromatic C-N coupling strategy. This is defined as a reaction occurring in the final steps of a synthesis (depth <= 3) that forms a C-N bond on an aromatic ring, where the product contains a specific heterocyclic scaffold. The function identifies the transformation by checking against a defined list of named reactions (e.g., Buchwald-Hartwig, Ullmann, SNAr) and verifies the product structure against a list of heterocycles of interest.
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

    # Track if we found the amination reaction
    found_amination = False
    # Track the depth at which amination occurs
    amination_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_amination, amination_depth, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata במקרה כזה, אם יש לך שאלות נוספות, אל תהסס לשאול."].get("rsmi", "")
            if not rsmi:
                return

            product_part = rsmi.split(">")[-1]

            # Check if this is a relevant C-N coupling reaction
            is_cn_coupling = False
            for reaction_name in AROMATIC_CN_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_cn_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            if is_cn_coupling:
                # Check if the product contains a heterocycle of interest
                has_heterocycle = False
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, product_part):
                        has_heterocycle = True
                        findings_json["atomic_checks"]["ring_systems"].append(ring)

                if has_heterocycle:
                    found_amination = True
                    # If multiple amination reactions are found, keep the one with the lowest depth
                    if amination_depth is None or depth < amination_depth:
                        amination_depth = depth
                    
                    # Record the co-occurrence constraint
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "scope": "single_reaction",
                            "targets": [
                                "AROMATIC_CN_COUPLING_REACTIONS",
                                "HETEROCYCLES_OF_INTEREST_IN_PRODUCT"
                            ]
                        }
                    })

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if amination was found at a late stage (depth <= 3)
    result = False
    if found_amination and amination_depth is not None and amination_depth <= 3:
        result = True
        # Record the positional constraint
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "The reaction satisfying the co-occurrence constraint",
                "position_type": "depth_from_root",
                "operator": "<=",
                "value": 3
            }
        })

    return result, findings_json
