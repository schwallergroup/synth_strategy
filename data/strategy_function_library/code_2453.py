from typing import Tuple, Dict, List
import copy
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
    This function detects if the synthetic route involves azide formation
    and subsequent transformations using the azide group.
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

    azide_molecules = set()  # Track molecules containing azide groups
    azide_reactions = []  # Track reactions involving azides
    azide_formation_reactions = []  # Track reactions where azide is formed
    
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal azide_reactions, azide_formation_reactions, azide_molecules, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains an azide group
            if checker.check_fg("Azide", mol_smiles):
                azide_molecules.add(mol_smiles)
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an azide formation reaction
            reactants_have_azide = any(checker.check_fg("Azide", r) for r in reactants)
            product_has_azide = checker.check_fg("Azide", product)

            # Azide formation: no azide in reactants, but azide in product
            if not reactants_have_azide and product_has_azide:
                azide_formation_reactions.append(depth)
                if "azide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("azide_formation")

            # Azide transformation: azide in reactants
            elif reactants_have_azide:
                azide_reactions.append(depth)
                if "azide_transformation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("azide_transformation")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if we have evidence of an azide-mediated strategy
    # Case 1: Azide is formed and then used in subsequent reactions
    if azide_formation_reactions and azide_reactions:
        result = True
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["azide_formation", "azide_transformation"]}})

    # Case 2: Starting material contains azide and it's used in multiple reactions
    if len(azide_molecules) > 0 and len(azide_reactions) >= 2:
        result = True
        # This constraint combines two aspects: presence of Azide FG and count of azide_transformation
        # We'll add both relevant constraints if they are not already present
        co_occurrence_azide_transform = {"type": "co-occurrence", "details": {"targets": ["Azide", "azide_transformation"]}}
        if co_occurrence_azide_transform not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(co_occurrence_azide_transform)
        
        count_azide_transform = {"type": "count", "details": {"target": "azide_transformation", "operator": ">=", "value": 2}}
        if count_azide_transform not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(count_azide_transform)

    # Case 3: Multiple azide transformations are present (even without explicit formation)
    if len(azide_reactions) >= 2:
        result = True
        count_azide_transform = {"type": "count", "details": {"target": "azide_transformation", "operator": ">=", "value": 2}}
        if count_azide_transform not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(count_azide_transform)

    return result, findings_json
