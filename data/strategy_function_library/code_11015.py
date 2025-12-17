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
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "N-methylation",
    "DMS Amine methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
]

SULFONAMIDE_FORMATION_REACTIONS = [
    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
    "Schotten-Baumann to ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a late-stage N-alkylation strategy where:
    1. The final reaction involves N-alkylation (particularly N-methylation)
    2. The route preserves sensitive functional groups like nitriles
    3. The synthesis follows a linear path with sulfonamide formation
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

    # Initialize tracking variables
    final_step_is_n_alkylation = False
    contains_sulfonamide_formation = False
    preserves_nitrile = True
    nitrile_present = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_n_alkylation, contains_sulfonamide_formation, preserves_nitrile, nitrile_present, findings_json

        if node["type"] == "mol":
            # Check if molecule contains nitrile
            if "smiles" in node:
                # Use checker function to detect nitrile
                if checker.check_fg("Nitrile", node["smiles"]):
                    nitrile_present = True
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    print(f"Found nitrile in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for nitrile preservation
            reactant_has_nitrile = False
            for reactant in reactants:
                if checker.check_fg("Nitrile", reactant):
                    reactant_has_nitrile = True
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    print(f"Found nitrile in reactant: {reactant}")
                    break

            product_has_nitrile = checker.check_fg("Nitrile", product)
            if product_has_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                print(f"Found nitrile in product: {product}")

            if reactant_has_nitrile and not product_has_nitrile:
                preserves_nitrile = False
                findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "destruction of Nitrile group"}})
                print("Nitrile not preserved in this reaction")

            # Check for N-alkylation in the final step (depth 1)
            if depth == 1:
                # Check for N-alkylation reactions using checker function
                for rxn_type in N_ALKYLATION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        final_step_is_n_alkylation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "N-alkylation", "position": "last_stage"}})
                        print(f"Final step is {rxn_type}")
                        break

            # Check for sulfonamide formation
            for rxn_type in SULFONAMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    contains_sulfonamide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Found sulfonamide formation: {rxn_type}")
                    break

            # Check if sulfonamide is formed in the reaction
            if not contains_sulfonamide_formation:
                reactants_have_sulfonamide = any(
                    checker.check_fg("Sulfonamide", r) for r in reactants
                )
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product)

                if not reactants_have_sulfonamide and product_has_sulfonamide:
                    contains_sulfonamide_formation = True
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    print("Found sulfonamide formation (functional group detected)")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'mol' node
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final step is N-alkylation: {final_step_is_n_alkylation}")
    print(f"Contains sulfonamide formation: {contains_sulfonamide_formation}")
    print(f"Preserves nitrile: {preserves_nitrile}")
    print(f"Nitrile present: {nitrile_present}")

    # The strategy is present if the final step is N-alkylation,
    # the route contains sulfonamide formation, and nitrile is preserved
    result = (
        final_step_is_n_alkylation
        and contains_sulfonamide_formation
        and preserves_nitrile
        and nitrile_present
    )

    if result:
        # Add the co-occurrence constraint if all conditions are met
        # This assumes that if the overall result is True, then the co-occurrence is also true.
        # The specific targets 'sulfonamide_formation' and 'Nitrile' are derived from the strategy JSON.
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["sulfonamide_formation", "Nitrile"]}})

    return result, findings_json
