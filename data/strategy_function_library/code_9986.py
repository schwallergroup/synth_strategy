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


N_ALKYLATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Reductive amination with aldehyde",
    "Reductive amination with ketone",
    "Reductive amination with alcohol",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "N-methylation",
    "Alkylation of amines",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Schotten-Baumann to ester",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid to amide conversion",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Ester with ammonia to amide",
    "Acylation of secondary amines with anhydrides",
    "{Schotten-Baumann_amide}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis pattern combining several key transformations. The pattern requires:
    1) at least one convergent step joining two complex fragments (defined as ring-containing);
    2) a nitro group reduction;
    3) an early-stage N-alkylation reaction from a defined list of reaction types;
    4) a late-stage amide formation reaction from a defined list of reaction types.
    All four conditions must be met somewhere in the synthetic route.
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

    # Initialize flags to track strategy components
    has_convergent_step = False
    has_nitro_reduction = False
    has_late_amide_formation = False
    has_early_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step, has_nitro_reduction, has_late_amide_formation, has_early_n_alkylation, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")

            # Check for convergent synthesis (multiple complex reactants)
            if len(reactants) >= 2:
                complex_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetRingInfo().NumRings() > 0:
                            complex_reactants += 1
                            if "any_ring" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("any_ring")
                    except:
                        pass

                if complex_reactants >= 2:
                    print(f"Found convergent step with complex reactants at depth {depth}")
                    has_convergent_step = True
                    if {"type": "count", "details": {"target": "ring-containing_reactants_in_one_step", "operator": ">=", "value": 2, "note": "At least one reaction in the route must have two or more reactants that contain a ring system."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "ring-containing_reactants_in_one_step", "operator": ">=", "value": 2, "note": "At least one reaction in the route must have two or more reactants that contain a ring system."}})

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro reduction at depth {depth}")
                has_nitro_reduction = True
                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                if {"type": "count", "details": {"target": "Reduction of nitro groups to amines", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "Reduction of nitro groups to amines", "operator": ">=", "value": 1}})

            # Check for early-stage N-alkylation (depth >= 3)
            if depth >= 3:  # Early stage
                for rxn in N_ALKYLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found early-stage N-alkylation at depth {depth}")
                        has_early_n_alkylation = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        if {"type": "positional", "details": {"target_group": "n_alkylation", "min_depth": 3, "note": "An N-alkylation reaction must occur at a depth of 3 or more."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target_group": "n_alkylation", "min_depth": 3, "note": "An N-alkylation reaction must occur at a depth of 3 or more."}})
                        break

            # Check for late-stage amide formation (depth <= 3)
            if depth <= 3:  # Late stage
                for rxn in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Found late-stage amide formation at depth {depth}")
                        has_late_amide_formation = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        if {"type": "positional", "details": {"target_group": "amide_formation", "max_depth": 3, "note": "An amide formation reaction must occur at a depth of 3 or less."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target_group": "amide_formation", "max_depth": 3, "note": "An amide formation reaction must occur at a depth of 3 or less."}})
                        break

        # Continue traversing children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Summary: Convergent step: {has_convergent_step}, Nitro reduction: {has_nitro_reduction}, "
        f"Late amide formation: {has_late_amide_formation}, Early N-alkylation: {has_early_n_alkylation}"
    )

    # Return True if all strategy components are present
    result = (
        has_convergent_step
        and has_nitro_reduction
        and has_late_amide_formation
        and has_early_n_alkylation
    )
    return result, findings_json
