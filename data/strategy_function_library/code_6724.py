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
    This function detects a strategy where an aromatic ether is formed by
    connecting an aromatic ring to a secondary carbon.
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

    has_aromatic_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_ether_formation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for relevant ether formation reactions
            is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)
            if is_williamson:
                if "Williamson Ether Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")

            is_mitsunobu = checker.check_reaction("Mitsunobu aryl ether", rsmi)
            if is_mitsunobu:
                if "Mitsunobu aryl ether" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu aryl ether")

            # Check if reactants contain a phenol
            has_aromatic_alcohol = False
            for reactant in reactants:
                if checker.check_fg("Phenol", reactant):
                    has_aromatic_alcohol = True
                    if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Phenol")
                    break

            # Check for aromatic ether in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Check for aromatic carbon connected to oxygen connected to any carbon
                aromatic_ether_pattern = Chem.MolFromSmarts("c-[O]-[C,c]")
                has_aromatic_ether = product_mol.HasSubstructMatch(aromatic_ether_pattern)

                # Check if the carbon connected to oxygen is secondary
                secondary_carbon_pattern = Chem.MolFromSmarts("c-O-[CH1](-[#6])(-[#6])")
                has_secondary_carbon = product_mol.HasSubstructMatch(secondary_carbon_pattern)

                # Verify the reaction is forming an aromatic ether with a secondary carbon
                # This logic primarily covers the Ar-OH + R-X pathway.
                if has_aromatic_alcohol and has_aromatic_ether:
                    if is_williamson or is_mitsunobu or has_secondary_carbon:
                        has_aromatic_ether_formation = True
                        # Record the 'aromatic_ether_formation' atomic check
                        if "aromatic_ether_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatic_ether_formation")
                        # Record the 'aromatic_ether_on_secondary_carbon' functional group if applicable
                        if has_secondary_carbon:
                            if "aromatic_ether_on_secondary_carbon" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("aromatic_ether_on_secondary_carbon")

                        # Record the structural constraint if all conditions are met
                        if has_aromatic_alcohol and has_aromatic_ether and (is_williamson or is_mitsunobu or has_secondary_carbon):
                            # This corresponds to the single structural constraint in the provided JSON
                            constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Phenol",
                                        "aromatic_ether_formation"
                                    ],
                                    "scope": "reaction_step",
                                    "condition": "A 'Phenol' must be present in the reactants and an aromatic ether must be formed in the product, and the reaction must be either a Williamson/Mitsunobu synthesis or result in a secondary carbon linkage."
                                }
                            }
                            if constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_aromatic_ether_formation, findings_json