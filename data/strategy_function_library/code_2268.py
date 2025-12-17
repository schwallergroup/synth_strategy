from typing import Tuple, Dict, List
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


KNOWN_MULTICOMPONENT_REACTIONS = [
    "Ugi reaction",
    "A3 coupling",
    "A3 coupling to imidazoles",
    "Petasis reaction with amines and boronic acids",
    "Petasis reaction with amines and boronic esters",
    "Petasis reaction with amines aldehydes and boronic acids",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for early-stage multicomponent reactions, defined as reactions occurring at any step except the last (depth >= 2) that either match a known named multicomponent reaction (including Ugi, A3, and Petasis variants) or involve three or more chemically distinct reactant molecules.
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

    found_multicomponent = False

    def dfs_traverse(node, depth=0):
        nonlocal found_multicomponent, findings_json

        if found_multicomponent:
            return  # Early return if we already found what we're looking for

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Check if this is an early-stage reaction (depth >= 2 means not the final step)
            if depth >= 2:
                # Record positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "multicomponent_reaction",
                        "position": "not_last_stage"
                    }
                })

                # Check for known multicomponent reactions first
                for name in KNOWN_MULTICOMPONENT_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        found_multicomponent = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        return

                # Count number of distinct reactants (at least 3 for multicomponent)
                if len(reactants) >= 3:
                    # Validate reactants are chemically distinct by converting to canonical SMILES
                    unique_reactants = set()
                    for r in reactants:
                        # Skip empty reactants
                        if not r.strip():
                            continue

                        try:
                            # RDKit handles atom-mapped SMILES directly. No manual cleaning needed.
                            mol = Chem.MolFromSmiles(r)
                            if mol:
                                canonical_smiles = Chem.MolToSmiles(mol)
                                unique_reactants.add(canonical_smiles)
                        except Exception:
                            # Skip invalid SMILES
                            continue

                    if len(unique_reactants) >= 3:
                        found_multicomponent = True
                        # Record count constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "distinct_reactants",
                                "operator": ">=",
                                "value": 3
                            }
                        })
                        return

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_multicomponent, findings_json
