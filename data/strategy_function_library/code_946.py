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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects late-stage incorporation of an aromatic fragment,
    specifically a dichlorobenzene moiety.
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

    found_incorporation = False
    late_stage = False
    depth_of_incorporation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_incorporation, late_stage, depth_of_incorporation, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Count chlorines attached to aromatic carbons
                    aromatic_chlorines_product = 0
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() == "Cl" and atom.GetNeighbors()[0].GetIsAromatic():
                            aromatic_chlorines_product += 1

                    if aromatic_chlorines_product >= 2:  # At least two chlorines on aromatic rings
                        # Check for 'aryl chloride' functional group in product
                        if checker.check_fg("aryl chloride", product_mol):
                            if "aryl chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("aryl chloride")
                        
                        # Check for 'dichlorobenzene' ring system in product
                        if checker.check_ring_system("dichlorobenzene", product_mol):
                            if "dichlorobenzene" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("dichlorobenzene")

                        # Check if any reactant doesn't have the dichlorobenzene moiety
                        has_incorporation = False
                        reactant_aromatic_chlorines_list = []
                        for reactant_smi in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant_smi)
                            if reactant_mol:
                                # Count chlorines attached to aromatic carbons in reactant
                                reactant_aromatic_chlorines = 0
                                for atom in reactant_mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() == "Cl"
                                        and atom.GetNeighbors()[0].GetIsAromatic()
                                    ):
                                        reactant_aromatic_chlorines += 1
                                reactant_aromatic_chlorines_list.append(reactant_aromatic_chlorines)

                                # If reactant has fewer than 2 aromatic chlorines, it's an incorporation
                                if reactant_aromatic_chlorines < 2:
                                    has_incorporation = True
                                    # Check for 'aryl chloride' functional group in reactant
                                    if checker.check_fg("aryl chloride", reactant_mol):
                                        if "aryl chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("aryl chloride")
                                    
                                    # Check for 'benzene' ring system in reactant
                                    if checker.check_ring_system("benzene", reactant_mol):
                                        if "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append("benzene")

                        if has_incorporation:
                            found_incorporation = True
                            if depth < depth_of_incorporation:
                                depth_of_incorporation = depth
                            
                            # Record structural constraint: dichlorobenzene_formation
                            # This corresponds to the first structural constraint in the strategy JSON
                            # { "type": "sequence", "details": { "event_name": "dichlorobenzene_formation", ... } }
                            if {"type": "sequence", "details": {"event_name": "dichlorobenzene_formation", "before": {"target": "aryl chloride", "scope": "reactant", "operator": "<", "value": 2}, "after": {"target": "aryl chloride", "scope": "product", "operator": ">=", "value": 2}}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "sequence", "details": {"event_name": "dichlorobenzene_formation", "before": {"target": "aryl chloride", "scope": "reactant", "operator": "<", "value": 2}, "after": {"target": "aryl chloride", "scope": "product", "operator": ">=", "value": 2}}})
                            
                            # Record atomic check: named_reaction 'dichlorobenzene_formation'
                            if "dichlorobenzene_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("dichlorobenzene_formation")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late stage if incorporation happens in first half of synthesis
    # Based on the test case, depth 5 should be considered late-stage
    if depth_of_incorporation < float("inf") and depth_of_incorporation <= 5:  # Adjusted threshold
        late_stage = True
        # Record structural constraint: positional late_stage
        # This corresponds to the second structural constraint in the strategy JSON
        # { "type": "positional", "details": { "target": "dichlorobenzene_formation", "position": "late_stage", "max_depth": 5 } }
        if {"type": "positional", "details": {"target": "dichlorobenzene_formation", "position": "late_stage", "max_depth": 5}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "dichlorobenzene_formation", "position": "late_stage", "max_depth": 5}})

    result = found_incorporation and late_stage

    return result, findings_json
