import copy
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
from typing import Tuple, Dict, List


def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage coupling reaction between two complex fragments, where 'complex' is defined as having more than 10 heavy atoms.
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

    complex_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal complex_coupling_detected, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Record late_stage constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "coupling_reaction",
                    "position": "late_stage"
                }
            })

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if there are at least 2 reactants
                if len(reactants) >= 2:
                    complex_fragments = 0

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if (
                                mol and mol.GetNumHeavyAtoms() > 10
                            ):  # Consider fragments with >10 heavy atoms as complex
                                complex_fragments += 1
                        except:
                            continue

                    if complex_fragments >= 2:
                        complex_coupling_detected = True
                        # Record complex_reactants constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants",
                                "operator": ">=",
                                "value": 2
                            }
                        })

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return complex_coupling_detected, findings_json
