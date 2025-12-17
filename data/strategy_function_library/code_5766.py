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


AMINE_ACYLATION_REACTION_TYPES = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with primary amine to imide",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of secondary amines with anhydrides",
    "Acylation of secondary amines",
    "Acylation of primary amines",
    "Schotten-Baumann_amide",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final step) amine acylation. The strategy is confirmed if the reaction matches a defined set of `AMINE_ACYLATION_REACTION_TYPES` or, alternatively, if a primary/secondary amine and an acylating agent are consumed to form an amide.
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

    final_step_is_acylation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_acylation, findings_json

        if final_step_is_acylation:
            return

        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_acylation_reaction = False
                for rxn_type in AMINE_ACYLATION_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_acylation_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if not is_acylation_reaction:
                    acylating_agents = ["Acyl halide", "Anhydride", "Carboxylic acid", "Ester"]
                    has_acylating_agent = False
                    for reactant in reactants_smiles:
                        for agent in acylating_agents:
                            if checker.check_fg(agent, reactant):
                                has_acylating_agent = True
                                if agent not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(agent)
                                break
                        if has_acylating_agent:
                            break

                    amine_types = ["Primary amine", "Secondary amine"]
                    has_amine_reactant = False
                    for reactant in reactants_smiles:
                        for amine_type in amine_types:
                            if checker.check_fg(amine_type, reactant):
                                has_amine_reactant = True
                                if amine_type not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(amine_type)
                                break
                        if has_amine_reactant:
                            break

                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                    has_amide_product = False
                    for amide_type in amide_types:
                        if checker.check_fg(amide_type, product_smiles):
                            has_amide_product = True
                            if amide_type not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(amide_type)
                            break

                    if has_acylating_agent and has_amine_reactant and has_amide_product:
                        is_acylation_reaction = True

                if is_acylation_reaction:
                    final_step_is_acylation = True
                    # Add the structural constraint if the final step is indeed an acylation
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "amine_acylation",
                            "position": "last_stage"
                        }
                    })

            except Exception:
                pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return final_step_is_acylation, findings_json
