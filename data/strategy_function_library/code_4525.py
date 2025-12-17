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
    Detects if the synthesis includes an intramolecular cyclization where a reactant containing a nitrile is converted into a product containing an amide and at least one new ring. This is a strong proxy for identifying nitrile-to-lactam formation.
    """
    found_nitrile_cyclization = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_cyclization, findings_json

        # For reaction nodes, check for nitrile to lactam conversion
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check reactants for nitriles
            reactants = reactants_part.split(".")
            has_nitrile_reactant = False
            for r in reactants:
                if checker.check_fg("Nitrile", r):
                    has_nitrile_reactant = True
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    break

            # Check product for amide
            product = Chem.MolFromSmiles(product_part)
            if product is None:
                for child in node.get("children", []):
                    # New logic: depth remains same if current node is 'reaction'
                    new_depth = depth if node["type"] == "reaction" else depth + 1
                    dfs_traverse(child, new_depth)
                return

            product_smiles = Chem.MolToSmiles(product)

            has_amide = False
            if checker.check_fg("Primary amide", product_smiles):
                has_amide = True
                if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
            if checker.check_fg("Secondary amide", product_smiles):
                has_amide = True
                if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
            if checker.check_fg("Tertiary amide", product_smiles):
                has_amide = True
                if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")

            product_rings = product.GetRingInfo().NumRings()

            # Use a check for intramolecularity as a key part of the strategy
            is_intramolecular = checker.check_reaction("Intramolecular amination", rsmi)
            if is_intramolecular:
                if "Intramolecular amination" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Intramolecular amination")

            # Count rings in reactants
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            reactant_rings = sum(
                mol.GetRingInfo().NumRings() for mol in reactant_mols if mol is not None
            )

            # Check for cyclization (more rings in product than in reactants)
            cyclization_occurred = product_rings > reactant_rings
            if cyclization_occurred:
                # This is a proxy for 'ring_formation' in the strategy
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Detect nitrile to lactam conversion
            if has_nitrile_reactant and has_amide and cyclization_occurred and is_intramolecular:
                found_nitrile_cyclization = True
                # Add the structural constraint if all conditions are met
                constraint_obj = {
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "Nitrile",
                            "Amide",
                            "ring_formation",
                            "Intramolecular amination"
                        ]
                    }
                }
                if constraint_obj not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint_obj)

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic: depth remains same if current node is 'reaction'
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_nitrile_cyclization, findings_json
