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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "{Schotten-Baumann_amide}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route uses a late-stage amide formation as the final step,
    specifically with an aromatic amine nucleophile. The reaction is identified by
    matching against a predefined list of amide formation reaction types.
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

    final_step_is_amide_formation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal final_step_is_amide_formation, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            is_final_step = current_depth == 1

            if is_final_step:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    reactant_smiles = [Chem.MolToSmiles(mol) for mol in reactant_mols if mol]

                    # Check for an aromatic amine (aniline) reactant
                    has_aromatic_amine = False
                    for r in reactant_smiles:
                        if checker.check_fg("Aniline", r):
                            has_aromatic_amine = True
                            if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Aniline")
                            break

                    # Check if this is a known amide formation reaction
                    is_amide_formation = False
                    for name in AMIDE_FORMATION_REACTIONS:
                        if checker.check_reaction(name, rsmi):
                            is_amide_formation = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            break

                    if is_amide_formation and has_aromatic_amine:
                        final_step_is_amide_formation = True
                        # Add structural constraints when the full condition is met
                        if {"type": "co-occurrence", "details": {"targets": ["amide_formation", "aniline_in_reactant"]}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amide_formation", "aniline_in_reactant"]}})
                        if {"type": "positional", "details": {"target": "amide_formation_with_aniline_reactant", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation_with_aniline_reactant", "position": "last_stage"}})

                except Exception as e:
                    # In a production environment, you might want to log this error.
                    # For this analysis, we'll let it pass silently to avoid crashing.
                    pass

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, current_depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)
    return final_step_is_amide_formation, findings_json
