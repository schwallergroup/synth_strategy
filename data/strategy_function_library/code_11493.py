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
    "pyridine",
    "pyrrole",
    "furan",
    "thiophene",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrazole",
    "triazole",
    "tetrazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "indole",
    "benzimidazole",
    "quinoline",
    "isoquinoline",
    "piperidine",
    "morpholine",
    "piperazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a heterocyclic scaffold is built and then elaborated through sequential functional group interconversions. The function first identifies the formation of a heterocycle from the predefined list `HETEROCYCLES_OF_INTEREST`. It then verifies that at least two sequential functional group interconversion (FGI) reactions occur after the heterocycle has been formed.
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

    # Track heterocycle presence and its depth
    heterocycle_found = False
    heterocycle_depth = -1
    found_heterocycle_name = None

    # Track functional group transformations
    fg_transformations = []

    def has_heterocycle(smiles):
        """Check if a molecule contains a heterocycle"""
        nonlocal findings_json, found_heterocycle_name
        try:
            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, smiles):
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    found_heterocycle_name = ring
                    print(f"Found heterocycle {ring} in {smiles}")
                    return True
            return False
        except Exception as e:
            print(f"Error checking heterocycle: {e}")
            return False

    def is_fg_transformation(rsmi):
        """Check if a reaction is a functional group transformation"""
        nonlocal findings_json
        try:
            # Check common functional group transformations
            fg_rxn_types = [
                "Esterification of Carboxylic Acids",
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Reduction of aldehydes and ketones to alcohols",
                "Oxidation of aldehydes to carboxylic acids",
                "Oxidation of alcohol to carboxylic acid",
                "Reduction of ester to primary alcohol",
                "Reduction of carboxylic acid to primary alcohol",
                "Nitrile to amide",
                "Reduction of nitrile to amine",
                "Reduction of primary amides to amines",
                "Reduction of secondary amides to amines",
                "Reduction of tertiary amides to amines",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Alcohol to chloride_SOCl2",
                "Alcohol to triflate conversion",
            ]

            for rxn_type in fg_rxn_types:
                if checker.check_reaction(rxn_type, rsmi):
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Found FG transformation: {rxn_type}")
                    return True

            return False
        except Exception as e:
            print(f"Error checking FG transformation: {e}")
            return False

    def dfs_traverse(node, depth=0, path=[]):
        nonlocal heterocycle_found, heterocycle_depth, fg_transformations, findings_json

        current_path = path + [node]

        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule contains a heterocycle
            if has_heterocycle(node["smiles"]):
                # Only count if not a starting material or if we haven't found a heterocycle yet
                if not node.get("in_stock", False) or not heterocycle_found:
                    heterocycle_found = True
                    heterocycle_depth = depth
                    print(f"Found heterocycle at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            reaction_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a functional group transformation
            if is_fg_transformation(reaction_smiles):
                # Only count transformations after the heterocycle is found
                if heterocycle_found:
                    fg_transformations.append((depth, reaction_smiles))
                    print(
                        f"Found functional group transformation at depth {depth}: {reaction_smiles}"
                    )

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth, current_path)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth to check if they're sequential
    fg_transformations.sort(key=lambda x: x[0])

    # Check if transformations are sequential after heterocycle formation
    sequential_transformations = False
    if len(fg_transformations) >= 2:
        # Only consider transformations after heterocycle is found
        post_heterocycle_transformations = [
            t for t in fg_transformations if t[0] >= heterocycle_depth
        ]

        if len(post_heterocycle_transformations) >= 2:
            depths = [d for d, _ in post_heterocycle_transformations]
            # Check if at least 2 transformations are close to each other
            sequential_count = 0
            for i in range(1, len(depths)):
                if depths[i] - depths[i - 1] <= 3:  # Allow for 3-step gap
                    sequential_count += 1

            sequential_transformations = (
                sequential_count >= 1
            )  # At least 2 sequential transformations

    print(f"Heterocycle found: {heterocycle_found}")
    print(f"Heterocycle depth: {heterocycle_depth}")
    print(f"Number of FG transformations: {len(fg_transformations)}")
    print(
        f"Post-heterocycle transformations: {len([t for t in fg_transformations if t[0] >= heterocycle_depth])}"
    )
    print(f"Sequential transformations: {sequential_transformations}")

    # Consider it a heterocycle elaboration if:
    # 1. We found a heterocycle
    # 2. At least 2 functional group transformations after heterocycle formation
    # 3. The transformations are sequential
    post_heterocycle_count = len([t for t in fg_transformations if t[0] >= heterocycle_depth])

    result = False
    if (
        heterocycle_found
        and post_heterocycle_count >= 2
        and sequential_transformations
    ):
        print("Detected heterocycle elaboration strategy")
        result = True
        # Add structural constraints to findings_json
        if heterocycle_found:
            # This implies 'ring_formation' occurred, though not explicitly a reaction name
            # We can represent this as the first event in the sequence
            if {"type": "sequence", "details": {"ordered_events": ["ring_formation", "functional_group_interconversion"]}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "sequence", "details": {"ordered_events": ["ring_formation", "functional_group_interconversion"]}})
        
        if post_heterocycle_count >= 2:
            if {"type": "count", "details": {"target": "functional_group_interconversion_after_ring_formation", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "functional_group_interconversion_after_ring_formation", "operator": ">=", "value": 2}})
        
        if sequential_transformations:
            if {"type": "count", "details": {"target": "sequential_fgi_pair_after_ring_formation", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "sequential_fgi_pair_after_ring_formation", "operator": ">=", "value": 1}})

    return result, findings_json
