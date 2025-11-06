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


BORYLATION_REACTIONS = [
    "Preparation of boronic acids",
    "Preparation of boronic acids without boronic ether",
    "Preparation of boronic acids from trifluoroborates",
    "Preparation of boronic ethers",
]

SUZUKI_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step strategy where a boronic acid or ester is formed and then consumed in a subsequent Suzuki coupling. The specific borylation and Suzuki reactions checked are enumerated in BORYLATION_REACTIONS and SUZUKI_REACTIONS.
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

    borylation_found = False
    borylation_depth = -1
    suzuki_coupling_found = False
    suzuki_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal borylation_found, borylation_depth, suzuki_coupling_found, suzuki_coupling_depth, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for borylation (preparation of boronic acids/esters)
                for rxn in BORYLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        borylation_found = True
                        borylation_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found borylation at depth {depth}, rsmi: {rsmi}")
                        break # Found one borylation reaction, no need to check others

                # Check for Suzuki coupling
                for rxn in SUZUKI_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        suzuki_coupling_found = True
                        suzuki_coupling_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        print(f"Found Suzuki coupling at depth {depth}, rsmi: {rsmi}")
                        break # Found one Suzuki reaction, no need to check others

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # In retrosynthetic analysis, borylation should have higher depth than Suzuki coupling
    # (borylation happens earlier in the forward synthesis, so appears deeper in retrosynthetic tree)
    strategy_present = (
        borylation_found
        and suzuki_coupling_found
        and borylation_depth
        > suzuki_coupling_depth  # In retrosynthesis, borylation appears at greater depth
    )

    if strategy_present:
        print(
            f"Found Suzuki coupling strategy with boronic acid intermediate: borylation at depth {borylation_depth}, coupling at depth {suzuki_coupling_depth}"
        )
        # Add the structural constraint if the strategy is found
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": {
                    "any_of": [
                        "Preparation of boronic acids",
                        "Preparation of boronic acids without boronic ether",
                        "Preparation of boronic acids from trifluoroborates",
                        "Preparation of boronic ethers"
                    ],
                    "description": "Borylation (preparation of boronic acid/ester)"
                },
                "after": {
                    "any_of": [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic acids OTf",
                        "Suzuki coupling with boronic esters",
                        "Suzuki coupling with boronic esters OTf",
                        "Suzuki coupling with sulfonic esters"
                    ],
                    "description": "Suzuki coupling"
                }
            }
        })
    else:
        if borylation_found and suzuki_coupling_found:
            print(
                f"Found both reactions but in wrong order: borylation at depth {borylation_depth}, coupling at depth {suzuki_coupling_depth}"
            )
        elif borylation_found:
            print(f"Found only borylation at depth {borylation_depth}")
        elif suzuki_coupling_found:
            print(f"Found only Suzuki coupling at depth {suzuki_coupling_depth}")
        else:
            print("Neither borylation nor Suzuki coupling found")

    return strategy_present, findings_json
