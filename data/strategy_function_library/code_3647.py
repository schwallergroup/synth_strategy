from typing import Tuple, Dict, List
import copy
from rdkit import Chem

# Refactoring for Enumeration: List moved to module-level constant
AMIDE_FORMATION_REACTION_TYPES = [
    "Carboxylic acid with primary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Nucleophiles by Carboxylic Acids",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects if the synthesis uses an early-stage amide formation, specifically looking for carboxylic acid + amine reaction."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    found_amide_formation = False
    target_molecule = None
    amide_forming_reactions = []

    # First pass to get the target molecule (depth 0)
    if route["type"] == "mol":
        target_molecule = Chem.MolFromSmiles(route["smiles"])
        print(f"Target molecule: {route['smiles']}")

    def check_which_reaction_types(rsmi):
        """Helper function to check which reaction types match"""
        matching_types = []
        # Uses the module-level constant
        for rxn_type in AMIDE_FORMATION_REACTION_TYPES:
            if checker.check_reaction(rxn_type, rsmi):
                matching_types.append(rxn_type)
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
        return matching_types

    # REMOVED redundant helper function check_amide_formation_pattern

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation, amide_forming_reactions, findings_json

        if node["type"] == "reaction" and depth >= 2:  # Early-stage reaction (depth 2 or higher)
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                if not rsmi:
                    return

                try:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product = Chem.MolFromSmiles(product_part)

                    if not product or not all(reactants):
                        print("Failed to parse reactants or product")
                        return

                    # Convert back to SMILES for checker functions
                    reactant_smiles = [Chem.MolToSmiles(mol) for mol in reactants]
                    product_smiles = Chem.MolToSmiles(product)

                    # Check for carboxylic acid in reactants
                    has_acid = any(
                        checker.check_fg("Carboxylic acid", smiles) for smiles in reactant_smiles
                    )
                    if has_acid and "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

                    # Check for ester in reactants (alternative to carboxylic acid)
                    has_ester = any(checker.check_fg("Ester", smiles) for smiles in reactant_smiles)
                    if has_ester and "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")

                    # Check for acyl halide in reactants (alternative to carboxylic acid)
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", smiles) for smiles in reactant_smiles
                    )
                    if has_acyl_halide and "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")

                    # Check for primary or secondary amine in reactants
                    has_primary_amine = any(checker.check_fg("Primary amine", smiles) for smiles in reactant_smiles)
                    has_secondary_amine = any(checker.check_fg("Secondary amine", smiles) for smiles in reactant_smiles)
                    has_amine = has_primary_amine or has_secondary_amine

                    if has_primary_amine and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if has_secondary_amine and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                    # Check for amide in product
                    has_primary_amide = checker.check_fg("Primary amide", product_smiles)
                    has_secondary_amide = checker.check_fg("Secondary amide", product_smiles)
                    has_tertiary_amide = checker.check_fg("Tertiary amide", product_smiles)
                    has_amide = has_primary_amide or has_secondary_amide or has_tertiary_amide

                    if has_primary_amide and "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                    if has_secondary_amide and "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                    if has_tertiary_amide and "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

                    # Check if this is specifically an amide formation reaction
                    matching_reaction_types = check_which_reaction_types(rsmi)
                    is_amide_formation = len(matching_reaction_types) > 0

                    # REMOVED redundant manual_check logic block

                    # MODIFIED conditional to fix logical flaw and correctly use both checks
                    if is_amide_formation or (
                        has_amide and ((has_acid or has_ester or has_acyl_halide) and has_amine)
                    ):
                        print(f"Found early-stage amide formation at depth {depth}")
                        amide_forming_reactions.append((product_smiles, depth))
                        found_amide_formation = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = False
    # If we found amide formation reactions, verify they're early stage
    if amide_forming_reactions:
        print(f"Found {len(amide_forming_reactions)} amide formation reactions")
        # Sort by depth (highest depth = earliest stage)
        amide_forming_reactions.sort(key=lambda x: x[1], reverse=True)
        print(f"Earliest amide formation at depth {amide_forming_reactions[0][1]}")

        # If the earliest amide formation is at depth 2 or higher, it's early stage
        if amide_forming_reactions[0][1] >= 2:
            result = True
            # Add the structural constraint if the condition is met
            if {"type": "positional", "details": {"target": "amide_formation", "position": "early_stage", "definition": "depth >= 2"}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "early_stage", "definition": "depth >= 2"}})

    return result, findings_json
