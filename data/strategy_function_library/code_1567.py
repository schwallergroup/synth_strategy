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


SUZUKI_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "{Suzuki}",
    "Suzuki",
]

BORYLATION_REACTIONS = [
    "Preparation of boronic acids",
    "Preparation of boronic ethers",
    "Preparation of boronic acids from trifluoroborates",
    "Preparation of boronic acids without boronic ether",
    "Synthesis of boronic acids",
]

AROMATIC_SUBSTITUTION_REACTIONS = [
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "Ullmann condensation",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation_heterocycles",
    "{N-arylation_heterocycles}",
    "{Buchwald-Hartwig}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects routes containing specific convergent synthesis patterns. The function identifies two main scenarios: (1) A late-stage Suzuki coupling (occurring within the last three steps), or (2) A combination of a borylation reaction with a subsequent aromatic substitution or nitro group reduction. The check for Suzuki couplings uses the following reaction names: ['Suzuki coupling with boronic acids', 'Suzuki coupling with boronic esters', 'Suzuki coupling with boronic acids OTf', 'Suzuki coupling with boronic esters OTf', 'Suzuki coupling with sulfonic esters', '{Suzuki}', 'Suzuki']. The check for borylation uses: ['Preparation of boronic acids', 'Preparation of boronic ethers', 'Preparation of boronic acids from trifluoroborates', 'Preparation of boronic acids without boronic ether', 'Synthesis of boronic acids']. The check for aromatic substitution includes SNAr, Ullmann, and Buchwald-Hartwig type reactions: ['heteroaromatic_nuc_sub', 'nucl_sub_aromatic_ortho_nitro', 'nucl_sub_aromatic_para_nitro', 'Ullmann condensation', 'Ullmann-Goldberg Substitution amine', 'Ullmann-Goldberg Substitution thiol', 'Ullmann-Goldberg Substitution aryl alcohol', 'N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)', 'Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine', 'Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine', 'N-arylation_heterocycles', '{N-arylation_heterocycles}', '{Buchwald-Hartwig}'].
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

    suzuki_coupling_found = False
    borylation_found = False
    snar_count = 0
    nitro_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, borylation_found, snar_count, nitro_reduction_found, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for Suzuki coupling (depth 0, 1, or 2 = late stage)
                if depth <= 2:
                    for r in SUZUKI_REACTIONS:
                        if checker.check_reaction(r, rsmi):
                            suzuki_coupling_found = True
                            if r not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r)
                            print(f"Found late-stage Suzuki coupling at depth {depth}")
                            break # Found one Suzuki reaction, no need to check others

                # Check for borylation
                for r in BORYLATION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        borylation_found = True
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        print("Found borylation reaction")
                        break # Found one borylation reaction, no need to check others

                # Check for SNAr, Ullmann, or Buchwald-Hartwig type reactions
                for r in AROMATIC_SUBSTITUTION_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        snar_count += 1
                        if r not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r)
                        print(f"Found aromatic substitution reaction (count: {snar_count})")
                        break # Found one aromatic substitution reaction, no need to check others

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_reduction_found = True
                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    print("Found nitro reduction")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'chemical', depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = suzuki_coupling_found or (
        borylation_found and (snar_count >= 1 or nitro_reduction_found)
    )

    # Add structural constraints if met
    if suzuki_coupling_found:
        # This corresponds to the positional constraint for Suzuki coupling
        suzuki_constraint = {
            "type": "positional",
            "details": {
                "target": "Suzuki coupling",
                "position": "last_three_stages"
            }
        }
        if suzuki_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(suzuki_constraint)

    if borylation_found and snar_count >= 1:
        # This corresponds to the co-occurrence constraint for borylation and aromatic substitution
        borylation_aromatic_constraint = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "borylation",
                    "aromatic_substitution"
                ]
            }
        }
        if borylation_aromatic_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(borylation_aromatic_constraint)

    if borylation_found and nitro_reduction_found:
        # This corresponds to the co-occurrence constraint for borylation and nitro reduction
        borylation_nitro_constraint = {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "borylation",
                    "Reduction of nitro groups to amines"
                ]
            }
        }
        if borylation_nitro_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(borylation_nitro_constraint)

    print(
        f"Strategy components found: Suzuki={suzuki_coupling_found}, Borylation={borylation_found}, SNAr count={snar_count}, Nitro reduction={nitro_reduction_found}"
    )
    print(f"Strategy present: {strategy_present}")

    return strategy_present, findings_json
