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


ETHER_FORMATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "Mitsunobu aryl ether",
    "Mitsunobu aryl ether (intramolecular)",
    "Alcohol to ether",
    "Chan-Lam etherification",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (depth <= 1) ether formation. This is triggered if a reaction results in a net increase of ether functional groups or if the reaction is one of the specified ETHER_FORMATION_REACTIONS, such as Williamson Ether Synthesis, Mitsunobu reactions, or Ullmann condensation.
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

    late_stage_ether = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ether, findings_json

        if node["type"] == "reaction" and depth <= 1:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if checker.check_fg("Ether", product):
                    if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")

                    # Check for known ether formation reaction types
                    is_ether_formation = False
                    for rxn in ETHER_FORMATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_ether_formation = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    # Count ethers in reactants
                    reactant_ether_count = 0
                    for reactant in reactants:
                        reactant_ether_count += sum(
                            1 for _ in checker.get_fg_atom_indices("Ether", reactant)
                        )

                    # Count ethers in product
                    product_ether_matches = sum(
                        1 for _ in checker.get_fg_atom_indices("Ether", product)
                    )

                    if product_ether_matches > reactant_ether_count:
                        if "net_ether_increase" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("net_ether_increase")

                    if product_ether_matches > reactant_ether_count or is_ether_formation:
                        late_stage_ether = True

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if late_stage_ether:
        # Add the structural constraint if late-stage ether formation is detected
        # This corresponds to the 'positional' constraint in the strategy JSON
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ether_formation",
                "position": "late_stage"
            }
        })

    return late_stage_ether, findings_json
