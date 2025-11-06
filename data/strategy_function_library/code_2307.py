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


BENZIMIDAZOLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_aldehyde",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_formation from aldehyde",
    "benzimidazole_formation from acyl halide",
    "benzimidazole_formation from ester/carboxylic acid",
    "benzimidazole_aldehyde",
]

SNAR_REACTIONS = [
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "heteroaromatic_nuc_sub",
    "N-arylation",
    "Ullmann-Goldberg Substitution amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Buchwald-Hartwig",
]

SUZUKI_REACTIONS = [
    "Suzuki",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a benzimidazole ring is formed
    in the final step of the synthesis, preceded by nitro reduction and early biaryl formation.
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

    # Initialize flags for key features
    benzimidazole_formation_at_depth_1 = False
    nitro_reduction = False
    suzuki_coupling = False
    snar_reaction = False

    # Track molecules through the synthesis to ensure continuity
    molecule_track = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal benzimidazole_formation_at_depth_1, nitro_reduction, suzuki_coupling, snar_reaction, max_depth, findings_json

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Store the relationship between product and reactants
            if product_smiles not in molecule_track:
                molecule_track[product_smiles] = {"depth": depth, "reactants": reactants_smiles}

            if depth == 1:
                # Check for benzimidazole formation in the final step
                if checker.check_ring("benzimidazole", product_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append("benzimidazole")
                    reactants_have_benzimidazole = any(
                        checker.check_ring("benzimidazole", r) for r in reactants_smiles
                    )
                    is_benzimidazole_formation = any(
                        checker.check_reaction(rxn, rsmi) for rxn in BENZIMIDAZOLE_FORMATION_REACTIONS
                    )
                    if not reactants_have_benzimidazole and is_benzimidazole_formation:
                        benzimidazole_formation_at_depth_1 = True
                        for rxn in BENZIMIDAZOLE_FORMATION_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        if {"type": "positional", "details": {"target": "ring_formation", "target_name": "benzimidazole", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "target_name": "benzimidazole", "position": "last_stage"}})

            else:
                # For reactions at other depths, check if the product is a reactant in a reaction at depth-1
                # This ensures continuity in the synthetic pathway
                is_part_of_pathway = False
                for mol, info in molecule_track.items():
                    if info["depth"] == depth - 1 and product_smiles in info["reactants"]:
                        is_part_of_pathway = True
                        break

                if is_part_of_pathway or depth == max_depth:
                    if depth == 3:
                        # Check for nitro reduction
                        if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                            nitro_reduction = True
                            if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                            if {"type": "positional", "details": {"target": "Reduction of nitro groups to amines", "position": {"operator": "==", "value": 3}}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Reduction of nitro groups to amines", "position": {"operator": "==", "value": 3}}})

                    elif depth == 5:
                        # Check for SNAr reaction
                        snar_found = False
                        for rxn in SNAR_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                snar_reaction = True
                                snar_found = True
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        if snar_found:
                            # Add the SNAR structural constraint if any SNAR reaction is found at depth 5
                            snar_constraint = {"type": "alternative", "details": {"description": "Requires one of the specified biaryl formation strategies to occur at the specified depth.", "options": [{"target_group": "SUZUKI_REACTIONS", "position": {"operator": ">=", "value": 3}}, {"target_group": "SNAR_REACTIONS", "position": {"operator": "==", "value": 5}}]}}
                            if snar_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(snar_constraint)

                    elif depth >= 3:
                        # Check for Suzuki coupling
                        suzuki_found = False
                        for rxn in SUZUKI_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                suzuki_coupling = True
                                suzuki_found = True
                                if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        if suzuki_found:
                            # Add the Suzuki structural constraint if any Suzuki reaction is found at depth >= 3
                            suzuki_constraint = {"type": "alternative", "details": {"description": "Requires one of the specified biaryl formation strategies to occur at the specified depth.", "options": [{"target_group": "SUZUKI_REACTIONS", "position": {"operator": ">=", "value": 3}}, {"target_group": "SNAR_REACTIONS", "position": {"operator": "==", "value": 5}}]}}
                            if suzuki_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(suzuki_constraint)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        benzimidazole_formation_at_depth_1
        and nitro_reduction
        and (suzuki_coupling or snar_reaction)  # Either Suzuki or SNAr should be present
    )

    return strategy_present, findings_json
