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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis that involves N-functionalization of a primary amine, where the amine is prepared via a route that includes either a Henry (nitroaldol) reaction or a separate nitro group reduction. This check identifies the presence of these key transformations anywhere within the synthetic tree.
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

    # Track if we found the key reactions
    found_nitroaldol = False
    found_nitro_reduction = False
    found_n_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitroaldol, found_nitro_reduction, found_n_functionalization, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitroaldol condensation (Henry reaction)
            if checker.check_reaction("Henry Reaction", rsmi):
                found_nitroaldol = True
                findings_json["atomic_checks"]["named_reactions"].append("Henry Reaction")
                # print(f"Found nitroaldol condensation at depth {depth}: {rsmi}")

            # Check for nitro reduction to amine
            nitro_group_in_reactants = False
            for r in reactants:
                if checker.check_fg("Nitro group", r):
                    nitro_group_in_reactants = True
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

            if nitro_group_in_reactants and not checker.check_fg("Nitro group", product):
                if checker.check_fg("Primary amine", product):
                    found_nitro_reduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("nitro_group_reduction")
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    # print(f"Found nitro reduction at depth {depth}: {rsmi}")

            # Check for N-functionalization of the amine
            primary_amine_in_reactants = False
            for r in reactants:
                if checker.check_fg("Primary amine", r):
                    primary_amine_in_reactants = True
                    if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

            if primary_amine_in_reactants:
                # Check for various N-functionalization products
                functionalization_fgs = [
                    "Primary amide", "Secondary amide", "Tertiary amide",
                    "Urea", "Sulfonamide", "Secondary amine", "Tertiary amine"
                ]
                for fg_name in functionalization_fgs:
                    if checker.check_fg(fg_name, product):
                        found_n_functionalization = True
                        findings_json["atomic_checks"]["named_reactions"].append("N-functionalization")
                        if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                        # print(f"Found N-functionalization at depth {depth}: {rsmi}")
                        break # Found one functionalization, no need to check others for this reaction

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary of findings
    # print(f"Nitroaldol condensation found: {found_nitroaldol}")
    # print(f"Nitro reduction found: {found_nitro_reduction}")
    # print(f"N-functionalization found: {found_n_functionalization}")

    # Determine overall result
    result = (found_nitroaldol or found_nitro_reduction) and found_n_functionalization

    # Add structural constraints if the overall result is True
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "The route must contain an N-functionalization step, and also either a Henry Reaction or a nitro group reduction step.",
                "targets": [
                    "N-functionalization",
                    [
                        "Henry Reaction",
                        "nitro_group_reduction"
                    ]
                ]
            }
        })

    # Return True if we found the key components of this strategy
    # Either nitroaldol OR nitro_reduction is required, along with N-functionalization
    return result, findings_json
