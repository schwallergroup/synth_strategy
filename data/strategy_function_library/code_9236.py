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


from rdkit import Chem

# Pre-compiled patterns for efficiency, defined at module level.
HETEROCYCLES_OF_INTEREST_SMARTS = [
    "[#6]1[#6][#8][#6][#7]1",  # oxazole
    "[#6]1[#6][#8][#6][#6]1",  # furan
    "[#6]1[#6][#7][#7][#6]1",  # pyrazole
]
HETEROCYCLE_PATTERNS = [Chem.MolFromSmarts(s) for s in HETEROCYCLES_OF_INTEREST_SMARTS]
HETEROCYCLE_NAMES = {
    "[#6]1[#6][#8][#6][#7]1": "oxazole",
    "[#6]1[#6][#8][#6][#6]1": "furan",
    "[#6]1[#6][#7][#7][#6]1": "pyrazole",
}

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects convergent synthesis involving heterocyclic fragments
    (oxazole, furan, pyrazole).
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

    convergent_reactions = 0
    heterocycle_involvement = False
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_reactions, heterocycle_involvement, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if we have multiple reactants (potential convergent synthesis)
                if len(reactants_smiles) > 1:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                    if all(r is not None for r in reactants):
                        # Check if any reactant has heterocycles
                        heterocycle_in_reactants = False
                        for r in reactants:
                            for i, p in enumerate(HETEROCYCLE_PATTERNS):
                                if r.HasSubstructMatch(p):
                                    heterocycle_in_reactants = True
                                    heterocycle_name = HETEROCYCLE_NAMES.get(HETEROCYCLES_OF_INTEREST_SMARTS[i])
                                    if heterocycle_name and heterocycle_name not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle_name)
                                    break
                            if heterocycle_in_reactants:
                                break

                        if heterocycle_in_reactants:
                            heterocycle_involvement = True
                            convergent_reactions += 1
                            if "convergent_step" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("convergent_step")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found convergent synthesis with heterocycles
    result = convergent_reactions >= 1 and heterocycle_involvement

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_reaction_with_heterocycle",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json
