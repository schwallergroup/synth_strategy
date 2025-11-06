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
    Detects a synthetic strategy involving a late-stage hydroxamic acid formation,
    where an alpha-beta unsaturated carbonyl system is present somewhere in the synthetic route.
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
    found_hydroxamic_acid_formation = False
    found_alpha_beta_unsaturated = False

    # Track the molecules through the synthesis path
    target_carboxylic_acid = None
    hydroxamic_acid_depth = float("inf")

    def is_hydroxamic_acid(smiles):
        """Check if a molecule contains a hydroxamic acid group"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Hydroxamic acid patterns: R-C(=O)-N(R')-OH or R-C(=O)-N(OH)-R'
                patterns = ["C(=O)N[OH]", "C(=O)N(O)"]
                for pattern in patterns:
                    if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern.replace("R", "*"))):
                        return True
            return False
        except Exception:
            return False

    def is_alpha_beta_unsaturated(smiles):
        """Check if a molecule contains an alpha-beta unsaturated system"""
        try:
            # Additional check for alpha-beta unsaturated carbonyl compounds
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Alpha-beta unsaturated carbonyl pattern
                pattern = "C=C[$(C=O)]"  # any carbonyl adjacent to C=C
                if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                    return True
            return False
        except Exception:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal found_hydroxamic_acid_formation
        nonlocal found_alpha_beta_unsaturated, target_carboxylic_acid, hydroxamic_acid_depth, findings_json

        # Check for molecule nodes
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check for hydroxamic acid in molecule
            if not found_hydroxamic_acid_formation and is_hydroxamic_acid(mol_smiles):
                found_hydroxamic_acid_formation = True
                hydroxamic_acid_depth = min(hydroxamic_acid_depth, depth)
                if "hydroxamic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("hydroxamic acid")

            # Check for alpha-beta unsaturated systems
            if not found_alpha_beta_unsaturated and is_alpha_beta_unsaturated(mol_smiles):
                found_alpha_beta_unsaturated = True
                if "alpha-beta unsaturated carbonyl" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("alpha-beta unsaturated carbonyl")

        # Check for reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                reactants = reactants_str.split(".")
                product_str = rsmi.split(">")[-1]

                # Check for hydroxamic acid formation
                if checker.check_reaction("Hydroxamic Synthesis", rsmi):
                    found_hydroxamic_acid_formation = True
                    hydroxamic_acid_depth = min(hydroxamic_acid_depth, depth)
                    if "Hydroxamic Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Hydroxamic Synthesis")

                    # Store the reactant (carboxylic acid) for continuity check
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            target_carboxylic_acid = reactant
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            break

                # Check for ester hydrolysis
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    if "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                    # Check if product is a carboxylic acid
                    if checker.check_fg("Carboxylic acid", product_str):
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        # Check if any reactant is an ester with alpha-beta unsaturated system
                        for reactant in reactants:
                            if checker.check_fg("Ester", reactant) and is_alpha_beta_unsaturated(
                                reactant
                            ):
                                found_alpha_beta_unsaturated = True
                                if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Ester")
                                if "alpha-beta unsaturated carbonyl" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("alpha-beta unsaturated carbonyl")
                                break

                # Check all molecules in the reaction for alpha-beta unsaturated systems
                if not found_alpha_beta_unsaturated:
                    # Check product
                    if is_alpha_beta_unsaturated(product_str):
                        found_alpha_beta_unsaturated = True
                        if "alpha-beta unsaturated carbonyl" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("alpha-beta unsaturated carbonyl")

                    # Check reactants
                    if not found_alpha_beta_unsaturated:
                        for reactant in reactants:
                            if is_alpha_beta_unsaturated(reactant):
                                found_alpha_beta_unsaturated = True
                                if "alpha-beta unsaturated carbonyl" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("alpha-beta unsaturated carbonyl")
                                break

            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if hydroxamic acid formation is at the late stage (depth 0 or 1)
    late_stage_hydroxamic = hydroxamic_acid_depth <= 1

    # Populate structural constraints based on final flags
    if found_hydroxamic_acid_formation and found_alpha_beta_unsaturated:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "Hydroxamic Synthesis",
                    "alpha-beta unsaturated carbonyl"
                ]
            }
        })
    
    if found_hydroxamic_acid_formation and late_stage_hydroxamic:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Hydroxamic Synthesis",
                "position": "last_stage"
            }
        })

    # Return True if all key features were found and hydroxamic acid formation is late-stage
    result = (
        found_hydroxamic_acid_formation and late_stage_hydroxamic and found_alpha_beta_unsaturated
    )
    return result, findings_json
