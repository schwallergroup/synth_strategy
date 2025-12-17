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
rng_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**rng_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


ACYLATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Friedel-Crafts acylation",
    "Schotten-Baumann to ester",
    "Acylation of secondary amines with anhydrides",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage (final 3 steps) introduction of a trifluoroacetyl group. This strategy is identified when a trifluoro-containing acylating agent participates in a reaction from the `ACYLATION_REACTIONS` list, resulting in a product that now contains a trifluoroacetyl group. The check ensures the group is newly introduced and was not present on all reactants.
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

    has_trifluoroacetyl_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoroacetyl_introduction, findings_json

        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 1-3)
            if 0 < depth <= 3:
                # Record positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "trifluoroacetyl_introduction_reaction",
                        "position": "late_stage",
                        "description": "The reaction must occur within the final 3 steps of the synthesis (depth > 0 and depth <= 3)."
                    }
                })
                try:
                    is_acylation = False
                    for rxn in ACYLATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_acylation = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)

                    has_trifluoro_in_product = checker.check_fg("Trifluoro group", product)
                    if has_trifluoro_in_product:
                        findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                    has_carbonyl_in_product = (
                        checker.check_fg("Ketone", product)
                        or checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                        or checker.check_fg("Ester", product)
                    )
                    if checker.check_fg("Ketone", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if checker.check_fg("Primary amide", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if checker.check_fg("Secondary amide", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                    if checker.check_fg("Ester", product):
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                    trifluoro_in_reactants = [
                        checker.check_fg("Trifluoro group", r) for r in reactants
                    ]

                    acylating_agents = [
                        (
                            checker.check_fg("Acyl halide", r)
                            or checker.check_fg("Anhydride", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Carboxylic acid", r)
                        )
                        for r in reactants
                    ]
                    for r_idx, r in enumerate(reactants):
                        if checker.check_fg("Acyl halide", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                        if checker.check_fg("Anhydride", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Anhydride")
                        if checker.check_fg("Ester", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        if checker.check_fg("Carboxylic acid", r):
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    trifluoroacetyl_agents = [
                        trifluoro_in_reactants[i] and acylating_agents[i]
                        for i in range(len(reactants))
                    ]

                    if (
                        is_acylation
                        and has_trifluoro_in_product
                        and has_carbonyl_in_product
                        and any(trifluoroacetyl_agents)
                        and not all(trifluoro_in_reactants)
                    ):
                        has_trifluoroacetyl_introduction = True
                        # Record co-occurrence constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "reaction_step",
                                "targets": [
                                    "acylation_reaction_from_list",
                                    "product_has_trifluoro_group",
                                    "product_has_carbonyl_group",
                                    "reactant_is_trifluoro_acylating_agent"
                                ],
                                "description": "A single reaction step must be an acylation from a predefined list, the product must contain trifluoro and carbonyl groups, and at least one reactant must be a trifluoro-acylating agent (i.e., contains both a trifluoro group and an acylating functional group)."
                            }
                        })
                        # Record negation constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "scope": "reaction_step",
                                "target": "all reactants contain 'Trifluoro group'",
                                "description": "The reaction is only valid if not all reactants contain a trifluoro group, ensuring the group is newly introduced to the main scaffold."
                            }
                        })

                except Exception as e:
                    pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    # Remove duplicate entries from atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return has_trifluoroacetyl_introduction, findings_json
