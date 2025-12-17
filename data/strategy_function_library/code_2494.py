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


from rdkit import Chem

# Refactoring for Enumeration: The list of reactions is moved to a module-level
# constant and corrected to include only valid olefination reactions.
OLEFINATION_REACTIONS = [
    "Wittig",
    "Julia Olefination",
    "Horner-Wadsworth-Emmons",
]

def main(route) -> Tuple[bool, Dict]:
    """This function detects if a late-stage olefination is used to form a cyano-substituted alkene (C=C-C#N). The olefination is identified as one of the following named reactions: Wittig, Julia Olefination, or Horner-Wadsworth-Emmons."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    olefination_with_nitrile_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal olefination_with_nitrile_detected, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")

            # Check if this is a late-stage reaction (depth ≤ 2)
            if depth <= 2:
                is_olefination = False

                # Check if it's a known olefination reaction from the corrected list.
                for reaction_type in OLEFINATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_olefination = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # The flawed fallback check for phosphonates has been removed.

                if is_olefination:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if product contains a nitrile group
                        if checker.check_fg("Nitrile", product):
                            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                            product_mol = Chem.MolFromSmiles(product)

                            # The redundant check for an alkene in the product has been removed.

                            # Check if at least one reactant contains a nitrile group
                            nitrile_in_reactants = False
                            for reactant in reactants:
                                if checker.check_fg("Nitrile", reactant):
                                    nitrile_in_reactants = True
                                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                                    break

                            if nitrile_in_reactants:
                                # Check if the cyano-substituted alkene group is newly formed
                                cyanoethylidene_pattern = Chem.MolFromSmarts("C=CC#N")

                                # Check if cyano-substituted alkene exists in product
                                if product_mol and product_mol.HasSubstructMatch(cyanoethylidene_pattern):
                                    # This implies 'formation_of_cyano_substituted_alkene' is detected
                                    if "formation_of_cyano_substituted_alkene" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("formation_of_cyano_substituted_alkene")

                                    # Check if cyano-substituted alkene exists in any reactant
                                    cyanoethylidene_in_reactants = False
                                    for reactant in reactants:
                                        reactant_mol = Chem.MolFromSmiles(reactant)
                                        if reactant_mol and reactant_mol.HasSubstructMatch(
                                            cyanoethylidene_pattern
                                        ):
                                            cyanoethylidene_in_reactants = True
                                            break

                                    if not cyanoethylidene_in_reactants:
                                        olefination_with_nitrile_detected = True
                                        # Add structural constraints if all conditions are met
                                        findings_json["structural_constraints"].append({
                                            "type": "co-occurrence",
                                            "details": {
                                                "scope": "single_reaction_step",
                                                "targets": [
                                                    [
                                                        "Wittig",
                                                        "Julia Olefination",
                                                        "Horner-Wadsworth-Emmons"
                                                    ],
                                                    "Nitrile",
                                                    "formation_of_cyano_substituted_alkene"
                                                ],
                                                "comment": "A single reaction step must be one of the specified olefinations, involve a Nitrile group (in both reactant and product), and result in the formation of a cyano-substituted alkene."
                                            }
                                        })
                                        findings_json["structural_constraints"].append({
                                            "type": "positional",
                                            "details": {
                                                "target": "olefination_forming_cyano_alkene",
                                                "position": "late_stage",
                                                "max_depth": 2,
                                                "comment": "The reaction step matching the co-occurrence criteria must occur at a depth of 2 or less from the final product."
                                            }
                                        })
                                # The flawed 'else' block with the alternative pattern has been removed.
                    except Exception:
                        # Silently ignore errors from RDKit parsing or SMILES splitting.
                        pass

        # Traverse children with new depth logic
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return olefination_with_nitrile_detected, findings_json
