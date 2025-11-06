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


AMIDE_FORMATION_FROM_ACYL_HALIDE_RXNS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage amide formation from an acyl halide (or analog) where the acyl component is substituted with a trifluoromethyl group.
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

    late_stage_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_found, findings_json

        if (
            node["type"] == "reaction" and depth <= 1
        ):
            # Positional constraint check
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "amide_formation_from_acyl_halide",
                    "position": "late_stage (depth <= 1)",
                    "description": "The qualifying amide formation reaction must occur at depth 0 or 1 in the synthesis tree."
                }
            })

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                is_amide_formation = False
                for reaction_type in AMIDE_FORMATION_FROM_ACYL_HALIDE_RXNS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_formation = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if not is_amide_formation:
                    return

                trifluoro_in_reactants = False
                for reactant in reactants:
                    if checker.check_fg("Trifluoro group", reactant):
                        trifluoro_in_reactants = True
                        if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
                        break

                trifluoro_in_product = checker.check_fg("Trifluoro group", product)
                if trifluoro_in_product and "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                if is_amide_formation and trifluoro_in_reactants and trifluoro_in_product:
                    for reactant in reactants:
                        if checker.check_fg("Acyl halide", reactant) and checker.check_fg(
                            "Trifluoro group", reactant
                        ):
                            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                            if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")

                            late_stage_amide_found = True
                            # Co-occurrence constraint check
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "amide_formation_from_acyl_halide",
                                        "Acyl halide",
                                        "Trifluoro group"
                                    ],
                                    "description": "An amide formation from an acyl halide must occur, where one of the reactants contains both an Acyl halide and a Trifluoro group. The four specific reaction names are treated as a single logical group."
                                }
                            })
                            break

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return late_stage_amide_found, findings_json
