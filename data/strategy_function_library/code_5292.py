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


HETEROCYCLES_OF_INTEREST = [
    "triazole", "pyridine", "pyrimidine", "pyrazine", "thiadiazole",
    "oxadiazole", "pyrazole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "isoxazole", "isothiazole", "tetrazole", "indole",
    "benzimidazole", "benzoxazole", "benzothiazole",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with primary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of secondary amines",
    "Acylation of primary amines",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage amide formation involving a reactant that contains a specified heterocycle. This strategy is triggered if the reaction matches a predefined list of amide-forming named reactions from the `AMIDE_FORMATION_REACTIONS` list and a reactant contains a heterocycle from the `HETEROCYCLES_OF_INTEREST` list.
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

    # Track if we found a valid late-stage amide formation with heterocycle
    found_valid_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_valid_reaction, findings_json

        # Check reaction nodes at late stages (depth 0, 1, or 2)
        if node["type"] == "reaction" and depth <= 2:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            is_amide_formation = False
            # Check for specific, reliable amide formation reaction types
            for rxn_type in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_amide_formation = True
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

            if is_amide_formation:
                reactant_has_heterocycle = False
                # Check if any reactant contains one of the specified heterocycles
                for r in reactants:
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, r):
                            reactant_has_heterocycle = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                if reactant_has_heterocycle:
                    found_valid_reaction = True
                    # Add structural constraints if both conditions are met
                    # This corresponds to the 'co-occurrence' constraint
                    if {"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["amide_formation", "heterocycle_in_reactant"], "description": "A reaction must be an amide formation (from the AMIDE_FORMATION_REACTIONS list) and one of its reactants must contain a heterocycle (from the HETEROCYCLES_OF_INTEREST list)."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"scope": "reaction_step", "targets": ["amide_formation", "heterocycle_in_reactant"], "description": "A reaction must be an amide formation (from the AMIDE_FORMATION_REACTIONS list) and one of its reactants must contain a heterocycle (from the HETEROCYCLES_OF_INTEREST list)."}})
                    # This corresponds to the 'positional' constraint
                    if {"type": "positional", "details": {"target": "amide_formation_with_heterocycle", "position": "late_stage", "condition": "depth <= 2", "description": "The qualifying reaction must occur within the final three steps of the synthesis (depth 0, 1, or 2)."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation_with_heterocycle", "position": "late_stage", "condition": "depth <= 2", "description": "The qualifying reaction must occur within the final three steps of the synthesis (depth 0, 1, or 2)."}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_valid_reaction, findings_json
