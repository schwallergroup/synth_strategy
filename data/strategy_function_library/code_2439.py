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


NITROGEN_PROTECTING_GROUPS = ["Phthalimide", "Tosyl", "Nosyl", "Acetyl", "Trifluoroacetyl"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses that use at least two different types of nitrogen protecting groups.
    It identifies the use of Boc via specific reaction templates, and the use of Cbz,
    Phthalimide, Tosyl, Nosyl, Acetyl, and Trifluoroacetyl by detecting their
    addition or removal in a reaction step.
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

    protecting_groups_used = set()

    def dfs_traverse(node, depth=0):
        nonlocal protecting_groups_used, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection reactions
                boc_protection_reactions = [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine"
                ]
                for rxn_name in boc_protection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        protecting_groups_used.add("Boc")
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)

                # Check for Boc deprotection reactions
                boc_deprotection_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine"
                ]
                for rxn_name in boc_deprotection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        protecting_groups_used.add("Boc")
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)

                # Check for Cbz protection/deprotection
                reactants_have_cbz = any(checker.check_fg("Cbz", r) for r in reactants)
                product_has_cbz = checker.check_fg("Cbz", product)

                if not reactants_have_cbz and product_has_cbz:
                    protecting_groups_used.add("Cbz")
                    if "Cbz" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Cbz")

                if reactants_have_cbz and not product_has_cbz:
                    protecting_groups_used.add("Cbz")
                    if "Cbz" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Cbz")

                # Check for other protection/deprotection events
                for pg in NITROGEN_PROTECTING_GROUPS:
                    reactants_have_pg = any(checker.check_fg(pg, r) for r in reactants)
                    product_has_pg = checker.check_fg(pg, product)

                    if (not reactants_have_pg and product_has_pg) or (reactants_have_pg and not product_has_pg):
                        protecting_groups_used.add(pg)
                        if pg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(pg)

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = len(protecting_groups_used) >= 2

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_nitrogen_protecting_groups",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
