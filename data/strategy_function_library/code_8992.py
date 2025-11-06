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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis uses a convergent approach where a
    trifluoromethyl-containing fragment is coupled in a late stage.
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

    # Track if we've found a late-stage coupling with a CF3-containing fragment
    late_stage_cf3_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cf3_coupling, findings_json

        if (node["type"] == "reaction" and depth <= 1 and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if one of the reactants contains a CF3 group
            cf3_reactants = 0
            for reactant in reactants:
                if checker.has_functional_group(reactant, "trifluoromethyl"):
                    cf3_reactants += 1
                    if "trifluoromethyl" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")

            # If we have multiple reactants and at least one contains CF3, it's a convergent coupling
            if len(reactants) >= 2 and cf3_reactants >= 1:
                late_stage_cf3_coupling = True
                if "convergent_step" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("convergent_step")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if late_stage_cf3_coupling:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "convergent_step",
                    "trifluoromethyl"
                ]
            }
        })
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "convergent_coupling_with_trifluoromethyl",
                "position": "last_stage"
            }
        })

    return late_stage_cf3_coupling, findings_json
