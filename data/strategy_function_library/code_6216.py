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
    This function detects a synthetic strategy involving formation of an aromatic ether
    with a cycloalkylmethyl group, specifically looking for aryl-O-CH2-cycloalkyl patterns.
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
    has_cycloalkylmethyl_ether = False

    def dfs_traverse(node, depth, max_depth):
        nonlocal has_cycloalkylmethyl_ether, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an ether formation reaction
            is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)
            is_mitsunobu = checker.check_reaction("Mitsunobu aryl ether", rsmi)

            if is_williamson or is_mitsunobu:
                if is_williamson:
                    if "Williamson Ether Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                if is_mitsunobu:
                    if "Mitsunobu aryl ether" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu aryl ether")

                # Check if product contains cycloalkylmethyl ether
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # SMARTS for aryl-O-CH2-cycloalkyl pattern (any cycloalkyl ring)
                    cycloalkyl_patterns = [
                        Chem.MolFromSmarts("c-O-C-C1CC1"),
                        Chem.MolFromSmarts("c-O-C-C1CCC1"),
                        Chem.MolFromSmarts("c-O-C-C1CCCC1"),
                        Chem.MolFromSmarts("c-O-C-C1CCCCC1"),
                        Chem.MolFromSmarts("c-O-C-C1CCCCCC1"),
                        Chem.MolFromSmarts("c-O-C-C1CCCCCCC1"),
                    ]

                    product_has_pattern = False
                    for pattern in cycloalkyl_patterns:
                        if product_mol.HasSubstructMatch(pattern):
                            product_has_pattern = True
                            break

                    if product_has_pattern:
                        # Check reactants: need a phenol and a cycloalkylmethyl compound
                        has_phenol = False
                        has_cycloalkylmethyl = False

                        for reactant in reactants:
                            # Check for phenol
                            if checker.check_fg("Phenol", reactant):
                                has_phenol = True
                                if "Phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Phenol")

                            # Check for cycloalkylmethyl with leaving group
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                cycloalkyl_lg_patterns = [
                                    Chem.MolFromSmarts("C1CC1C[F,Cl,Br,I,O]"),
                                    Chem.MolFromSmarts("C1CCC1C[F,Cl,Br,I,O]"),
                                    Chem.MolFromSmarts("C1CCCC1C[F,Cl,Br,I,O]"),
                                    Chem.MolFromSmarts("C1CCCCC1C[F,Cl,Br,I,O]"),
                                    Chem.MolFromSmarts("C1CCCCCC1C[F,Cl,Br,I,O]"),
                                    Chem.MolFromSmarts("C1CCCCCCC1C[F,Cl,Br,I,O]"),
                                ]

                                for pattern in cycloalkyl_lg_patterns:
                                    if reactant_mol.HasSubstructMatch(pattern):
                                        has_cycloalkylmethyl = True
                                        break

                        if has_phenol and has_cycloalkylmethyl:
                            has_cycloalkylmethyl_ether = True
                            # Add structural constraint if all conditions are met
                            constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Williamson Ether Synthesis",
                                        "Mitsunobu aryl ether",
                                        "Phenol"
                                    ]
                                }
                            }
                            if constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(constraint_obj)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth, max_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1, max_depth)

    # Start traversal from the root, assuming max_depth is calculated by the environment.
    dfs_traverse(route, 1, 0) # The 0 is a placeholder for a calculated max_depth.

    return has_cycloalkylmethyl_ether, findings_json
