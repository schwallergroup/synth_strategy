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


# Refactored list of amidation reactions
AMIDATION_REACTIONS_OF_INTEREST = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects indole formation via azide reduction/cyclization followed by late-stage amidation.
    This strategy involves:
    1. Use of azide intermediate
    2. Indole formation
    3. Late-stage amidation
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

    has_azide = False
    has_indole_formation = False
    has_late_amidation = False

    # Track molecules with azides to ensure connection to indole formation
    azide_containing_molecules = set()
    indole_containing_molecules = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_azide, has_indole_formation, has_late_amidation, findings_json

        if node["type"] == "mol":
            # Check if this molecule contains an azide
            if checker.check_fg("Azide", node["smiles"]):
                has_azide = True
                azide_containing_molecules.add(node["smiles"])
                findings_json["atomic_checks"]["functional_groups"].append("Azide")
                print(f"Found azide at depth {depth}: {node['smiles']}")

            # Track indole-containing molecules
            if checker.check_ring("indole", node["smiles"]):
                indole_containing_molecules.add(node["smiles"])
                findings_json["atomic_checks"]["ring_systems"].append("indole")
                print(f"Found indole-containing molecule at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for indole formation
            reactants_have_indole = any(checker.check_ring("indole", r) for r in reactants_smiles)
            product_has_indole = checker.check_ring("indole", product_smiles)

            # Check if this reaction forms an indole from an azide
            if not reactants_have_indole and product_has_indole:
                # Check if any reactant contains an azide
                if any(checker.check_fg("Azide", r) for r in reactants_smiles) or any(
                    r in azide_containing_molecules for r in reactants_smiles
                ):
                    has_indole_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation") # Assuming 'ring_formation' covers indole formation
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "Azide", "after": "indole_formation", "description": "The indole ring must be formed from a reactant containing an azide group."}})
                    print(f"Found indole formation at depth {depth}: {rsmi}")

            # Check for amidation in the late stages (depth 0, 1, or 2)
            if depth <= 2:
                # Check for various amidation reactions
                amidation_found_in_reaction = False
                for rxn in AMIDATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn, rsmi):
                        amidation_found_in_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if amidation_found_in_reaction:
                    # Check if the product or any reactant contains an indole
                    # This ensures the amidation is related to the indole structure
                    if product_has_indole or reactants_have_indole:
                        has_late_amidation = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amidation", "position": "depth <= 2", "description": "Amidation must occur in the late stages of the synthesis (depth 0, 1, or 2)."}})
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amidation", "indole"], "scope": "single_reaction", "description": "The amidation reaction must involve a molecule that contains an indole ring."}})
                        print(f"Found late-stage amidation at depth {depth}: {rsmi}")

                # Additional check: if product contains both indole and amide but reactants don't have amide
                product_has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                product_has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                product_has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)

                if product_has_primary_amide: findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                if product_has_secondary_amide: findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                if product_has_tertiary_amide: findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                if product_has_indole and (
                    product_has_primary_amide
                    or product_has_secondary_amide
                    or product_has_tertiary_amide
                ):
                    reactants_have_amide = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants_smiles
                    )
                    if not reactants_have_amide:
                        has_late_amidation = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amidation", "position": "depth <= 2", "description": "Amidation must occur in the late stages of the synthesis (depth 0, 1, or 2)."}})
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["amidation", "indole"], "scope": "single_reaction", "description": "The amidation reaction must involve a molecule that contains an indole ring."}})
                        print(f"Found late-stage amide formation at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if all strategy components are found
    result = has_azide and has_indole_formation and has_late_amidation

    if result:
        # Add the overall co-occurrence constraint if all conditions are met
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Azide", "indole_formation", "amidation"], "description": "The overall strategy requires the presence of an azide intermediate, an indole formation step, and an amidation step."}})

    print(
        f"Strategy components: Azide: {has_azide}, Indole formation: {has_indole_formation}, Late amidation: {has_late_amidation}"
    )
    return result, findings_json
