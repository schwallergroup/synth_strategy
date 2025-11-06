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


CARBOXYLIC_ACID_PROTECTION_REACTIONS = [
    "Protection of carboxylic acid",
    "Esterification of Carboxylic Acids",
    "O-alkylation of carboxylic acids with diazo compounds",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage (depth <= 1) protection
    of a carboxylic acid. This is identified by checking for specific reaction
    types, including 'Protection of carboxylic acid', 'Esterification of
    Carboxylic Acids', and 'O-alkylation of carboxylic acids with diazo
    compounds', and confirming the consumption of a carboxylic acid functional group.
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

    protection_at_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal protection_at_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                has_acid_in_reactants = False
                for smi in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", smi):
                        has_acid_in_reactants = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                has_ester_in_product = False
                if checker.check_fg("Ester", product_smiles):
                    has_ester_in_product = True
                    if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                is_protection_reaction = False
                for name in CARBOXYLIC_ACID_PROTECTION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_protection_reaction = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

                if has_acid_in_reactants and has_ester_in_product and is_protection_reaction:
                    acid_count_reactants = sum(
                        len(checker.get_fg_atom_indices("Carboxylic acid", smi))
                        for smi in reactants_smiles
                    )
                    acid_count_product = len(
                        checker.get_fg_atom_indices("Carboxylic acid", product_smiles)
                    )

                    if acid_count_reactants > acid_count_product:
                        protection_at_depth = depth
                        # Record the co-occurrence constraint
                        if {"type": "co-occurrence", "details": {"description": "A reaction step must be a named protection reaction, have a 'Carboxylic acid' in the reactants, and an 'Ester' in the product.", "targets": ["Protection of carboxylic acid", "Carboxylic acid", "Ester"], "scope": "single_reaction_step"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "description": "A reaction step must be a named protection reaction, have a 'Carboxylic acid' in the reactants, and an 'Ester' in the product.",
                                    "targets": [
                                        "Protection of carboxylic acid",
                                        "Carboxylic acid",
                                        "Ester"
                                    ],
                                    "scope": "single_reaction_step"
                                }
                            })

            except Exception:
                # Silently ignore errors in reaction analysis to avoid crashing the whole run
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # If current node is 'reaction', depth remains the same for children (which are chemicals)
            # If current node is 'chemical', depth increases for children (which are reactions)
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it happens at depth 0 or 1
    is_late_stage = protection_at_depth is not None and protection_at_depth <= 1

    if is_late_stage:
        # Record the positional constraint
        if {"type": "positional", "details": {"target": "A named carboxylic acid protection reaction where a 'Carboxylic acid' functional group is consumed to form an 'Ester'.", "position": "late_stage (depth <= 1)"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "A named carboxylic acid protection reaction where a 'Carboxylic acid' functional group is consumed to form an 'Ester'.",
                    "position": "late_stage (depth <= 1)"
                }
            })

    return is_late_stage, findings_json
