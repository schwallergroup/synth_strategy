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
    This function detects a convergent coupling strategy where two complex fragments
    are combined in a late-stage reaction.
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

    late_stage_convergent_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_convergent_coupling, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late-stage reaction
            # Record late_stage positional constraint
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "any_reaction",
                    "position": "late_stage"
                }
            })

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if there are at least 2 complex fragments
                complex_fragments = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Define complexity as having at least 15 heavy atoms
                            if mol.GetNumHeavyAtoms() >= 15:
                                complex_fragments += 1
                    except:
                        continue

                if complex_fragments >= 2:
                    late_stage_convergent_coupling = True
                    # Record count constraint
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactant_with_heavy_atom_count_>=_15",
                            "operator": ">=",
                            "value": 2
                        }
                    })

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_convergent_coupling, findings_json