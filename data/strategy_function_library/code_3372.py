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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Schotten-Baumann_amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis featuring a late-stage amide coupling. The coupling is identified if it matches a specific list of reaction types, including various acylations and named reactions like Schotten-Baumann.
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

    # Track if we found a late-stage amide coupling of complex fragments
    found_late_stage_amide_coupling = False

    def is_complex_fragment(smiles):
        """Check if a fragment is complex (has rings or multiple functional groups)"""
        # A fragment is complex if it has at least one ring.
        if checker.get_ring_count(smiles) > 0:
            if "Ring system" not in findings_json["atomic_checks"]["ring_systems"]:
                # This is a generic check, so we'll just add a placeholder if a ring is found
                # The original strategy JSON doesn't list specific ring systems, so we can't add a specific one.
                # For this function, we'll just note that a ring system was detected as part of complexity.
                pass # No specific ring system name to add from the strategy JSON
            return True

        # Or if it has at least two functional groups.
        fg_count = 0
        fg_types = [
            "Carboxylic acid",
            "Ester",
            "Primary amide",
            "Secondary amide",
            "Tertiary amide",
            "Primary amine",
            "Secondary amine",
            "Tertiary amine",
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Aldehyde",
            "Ketone",
            "Nitrile",
            "Nitro group",
            "Ether",
            "Alkyne",
            "Alkene",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
            "Aromatic halide",
            "Azide",
            "Sulfonamide",
            "Boronic acid",
            "Boronic ester",
            "Phenol",
            "Thiol",
            "Sulfone",
            "Sulfoxide",
        ]

        for fg in fg_types:
            if checker.check_fg(fg, smiles):
                fg_count += 1
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                if fg_count >= 2:
                    return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_amide_coupling, findings_json

        # Process reaction nodes
        if node["type"] == "reaction" and node.get("children"):
            # Check if this is a late-stage reaction (depth 0 or 1)
            is_late_stage = (depth <= 1)
            if is_late_stage:
                if {"type": "positional", "details": {"target": "Amide Coupling", "position": "late_stage (depth <= 1)"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Amide Coupling", "position": "late_stage (depth <= 1)"}})

                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_smiles = rsmi.split(">")[0].split(".")

                    # Check if it's an amide coupling reaction by iterating through a defined list.
                    is_amide_coupling = False
                    for reaction_name in AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            is_amide_coupling = True
                            if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            break

                    # Check if complex fragments are being joined
                    if is_amide_coupling:
                        complex_reactants = 0
                        for r in reactants_smiles:
                            if is_complex_fragment(r):
                                complex_reactants += 1

                        # A convergent step joins at least two complex fragments.
                        is_convergent = complex_reactants >= 2

                        if is_convergent:
                            if {"type": "count", "details": {"target": "complex_reactants_in_amide_coupling", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactants_in_amide_coupling", "operator": ">=", "value": 2}})
                            found_late_stage_amide_coupling = True

                except (KeyError, Exception):
                    # Errors in reaction processing are ignored
                    pass

        # Continue traversing
        for child in node.get("children", []):
            # Stop traversing if we've already found the feature
            if found_late_stage_amide_coupling:
                break
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_late_stage_amide_coupling, findings_json
