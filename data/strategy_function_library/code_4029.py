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
    Detects if an O-demethylation reaction occurs within the last four steps of the synthesis. It first checks for a named reaction ('Cleavage of methoxy ethers to alcohols') and then falls back to a structural check for a decrease in methoxy groups and a corresponding increase in hydroxyl groups.
    """
    o_demethylation_detected = False

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
        nonlocal o_demethylation_detected, findings_json

        if node["type"] == "reaction" and 1 <= depth <= 4:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi):
                    print(f"O-demethylation reaction detected at depth {depth}")
                    o_demethylation_detected = True
                    findings_json["atomic_checks"]["named_reactions"].append("Cleavage of methoxy ethers to alcohols")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "O-demethylation", "position": "within_last_4_stages"}})
                    return

                reactant_has_ether = False
                decrease_in_methoxy_count = False
                increase_in_hydroxyl_count = False

                for reactant_smiles in reactants_smiles:
                    if checker.check_fg("Ether", reactant_smiles):
                        reactant_has_ether = True
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")

                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if reactant_mol and product_mol:
                            methoxy_pattern = Chem.MolFromSmarts("[CH3][OX2][#6]")
                            hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")

                            r_methoxy_matches = reactant_mol.GetSubstructMatches(methoxy_pattern)
                            p_hydroxyl_matches = product_mol.GetSubstructMatches(hydroxyl_pattern)

                            r_methoxy_count = len(r_methoxy_matches)
                            p_methoxy_count = len(product_mol.GetSubstructMatches(methoxy_pattern))
                            r_hydroxyl_count = len(
                                reactant_mol.GetSubstructMatches(hydroxyl_pattern)
                            )
                            p_hydroxyl_count = len(p_hydroxyl_matches)

                            if p_methoxy_count < r_methoxy_count:
                                decrease_in_methoxy_count = True
                                if "methoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("methoxy")

                            if p_hydroxyl_count > r_hydroxyl_count:
                                increase_in_hydroxyl_count = True
                                if "hydroxyl" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("hydroxyl")

                            if (
                                reactant_has_ether
                                and decrease_in_methoxy_count
                                and increase_in_hydroxyl_count
                            ):
                                print(
                                    f"O-demethylation functional group change detected at depth {depth}"
                                )
                                print(f"Reactant: {reactant_smiles}")
                                print(f"Product: {product_smiles}")
                                o_demethylation_detected = True
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "O-demethylation", "position": "within_last_4_stages"}})
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["reactant_has_ether", "decrease_in_methoxy_count", "increase_in_hydroxyl_count"], "scope": "reaction_step"}})
                                return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return o_demethylation_detected, findings_json
