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
    This function detects if the synthetic route involves protection of a carboxylic acid
    as a tert-butyl ester in the early stages of synthesis.
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

    protection_found = False

    def dfs_traverse(node, depth):
        nonlocal protection_found, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a carboxylic acid group
                has_carboxylic_acid_reactant = False
                for reactant in reactants:
                    if reactant and checker.check_fg("Carboxylic acid", reactant):
                        has_carboxylic_acid_reactant = True
                        if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        break

                # Check if product has a tert-butyl ester group
                has_tert_butyl_ester = False
                if product:
                    # Check specifically for tert-butyl group
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        tert_butyl_pattern = Chem.MolFromSmarts(
                            "[C;H0]([C;H3])([C;H3])([C;H3])[O;D2][C;D3](=[O;D1])"
                        )
                        if mol.HasSubstructMatch(tert_butyl_pattern):
                            has_tert_butyl_ester = True
                            if "tert-butyl ester" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("tert-butyl ester")

                # Check if this is a protection reaction
                is_protection_reaction = False
                if has_carboxylic_acid_reactant and has_tert_butyl_ester:
                    if checker.check_reaction(
                        "Esterification of Carboxylic Acids", rsmi
                    ):
                        is_protection_reaction = True
                        if "Esterification of Carboxylic Acids" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                    elif checker.check_reaction("Protection of carboxylic acid", rsmi):
                        is_protection_reaction = True
                        if "Protection of carboxylic acid" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Protection of carboxylic acid")

                # If we found a protection reaction at high depth (early in synthesis)
                # Early stage is typically at depth >= 2 (further from target)
                if is_protection_reaction and depth >= 2:
                    protection_found = True
                    print(f"Found tert-butyl ester protection at depth: {depth}")
                    # Add the structural constraint if detected
                    if {"type": "positional", "details": {"target": "Protection of carboxylic acid as tert-butyl ester", "position": "early_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Protection of carboxylic acid as tert-butyl ester", "position": "early_stage"}})

        # Continue traversing with increased depth
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal with depth 0
    dfs_traverse(route, 0)
    return protection_found, findings_json