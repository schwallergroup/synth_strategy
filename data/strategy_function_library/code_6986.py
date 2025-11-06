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
    Detects a strategy where a named heterocycle formation reaction occurs in the final or penultimate step,
    and a sulfonamide group is present as an intermediate or is formed at some point in the synthesis.
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

    # Track if we found the key features
    found_final_heterocycle_formation = False
    found_sulfonamide_intermediate = False

    # Define the list of heterocycle formation reactions for easier management
    HETEROCYCLE_FORMATION_REACTIONS = [
        "Formation of NOS Heterocycles",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{thiazole}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{pyrazole}",
        "{oxadiazole}",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "{Paal-Knorr pyrrole}",
        "Paal-Knorr pyrrole synthesis",
        "{Fischer indole}",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)"
    ]

    SULFONAMIDE_SYNTHESIS_REACTIONS = [
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine"
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_final_heterocycle_formation, found_sulfonamide_intermediate, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is the final or penultimate step (depth 0 or 1)
            if depth <= 1:
                print(f"Checking for heterocycle formation at depth {depth}")
                # Check for heterocycle formation reactions
                heterocycle_formation = False

                # Check common heterocycle formation reactions
                for reaction_name in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        heterocycle_formation = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        print(f"Found heterocycle formation reaction: {reaction_name} at depth {depth}")

                # Set the flag if we found heterocycle formation at depth 0 or 1
                if heterocycle_formation:
                    if depth == 0:
                        found_final_heterocycle_formation = True
                        print("Confirmed heterocycle formation in final step")
                    elif depth == 1 and not found_final_heterocycle_formation:
                        # Only set if not already found at depth 0
                        found_final_heterocycle_formation = True
                        print("Confirmed heterocycle formation in penultimate step")

            # Check for sulfonamide intermediate
            sulfonamide_found_in_reactants = False
            for r in reactants_smiles:
                if checker.check_fg("Sulfonamide", r):
                    sulfonamide_found_in_reactants = True
                    if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                    break

            if sulfonamide_found_in_reactants or checker.check_fg("Sulfonamide", product_smiles):
                found_sulfonamide_intermediate = True
                if "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")
                print(f"Found sulfonamide at depth {depth}")

            # Check for sulfonamide formation reaction
            for reaction_name in SULFONAMIDE_SYNTHESIS_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    found_sulfonamide_intermediate = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    print(f"Found sulfonamide formation reaction: {reaction_name} at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary of findings
    print("\nSummary of findings:")
    print(f"Heterocycle formation in final/penultimate step: {found_final_heterocycle_formation}")
    print(f"Sulfonamide intermediate: {found_sulfonamide_intermediate}")

    # Strategy is present if we found heterocycle formation in the final step
    # and the sulfonamide key feature
    result = found_final_heterocycle_formation and found_sulfonamide_intermediate

    # Add structural constraints if conditions are met
    if found_final_heterocycle_formation:
        # This corresponds to the 'positional' constraint for heterocycle formation
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target_type": "named_reaction",
                "targets": HETEROCYCLE_FORMATION_REACTIONS,
                "position": "last_or_penultimate_stage"
            }
        })

    if found_sulfonamide_intermediate:
        # This corresponds to the 'count' constraint for sulfonamide presence/formation
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target_type": "any",
                "targets": ["Sulfonamide"] + SULFONAMIDE_SYNTHESIS_REACTIONS,
                "operator": ">=",
                "value": 1
            }
        })

    print(f"Final result: {result}")
    return result, findings_json
