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
    Detects late-stage (final step) etherification of an indole core. The strategy is positive if an indole-containing reactant acts as either the oxygen nucleophile (e.g., hydroxy-indole) or the electrophile (e.g., halo-indole) in a C-O bond-forming reaction.
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

    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is the final reaction (depth 0 or 1)
            node_depth = node["metadata"].get("depth", depth)

            if node_depth <= 1:
                # Check for indole in reactants
                indole_in_reactants = False
                indole_reactant = None
                for reactant in reactants_smiles:
                    try:
                        if checker.check_ring("indole", reactant):
                            indole_in_reactants = True
                            indole_reactant = reactant
                            if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("indole")
                            break
                    except Exception:
                        pass

                # Check if product contains indole
                indole_in_product = False
                try:
                    indole_in_product = checker.check_ring("indole", product_smiles)
                    if indole_in_product and "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("indole")
                except Exception:
                    pass

                if indole_in_reactants and indole_reactant and indole_in_product:
                    # Structural constraint: indole and C-O bond formation co-occurrence
                    if {"type": "co-occurrence", "details": {"targets": ["indole", "C-O bond formation"]}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["indole", "C-O bond formation"]}})

                    # Check if this is a C-O bond forming reaction
                    co_bond_reaction = False

                    # Check for common C-O bond forming reactions
                    co_bond_reactions_list = [
                        "Williamson Ether Synthesis",
                        "Mitsunobu aryl ether",
                        "Chan-Lam etherification",
                        "Ullmann-Goldberg Substitution aryl alcohol",
                        "Alcohol to ether",
                        "nucl_sub_aromatic_ortho_nitro",
                        "nucl_sub_aromatic_para_nitro",
                        "heteroaromatic_nuc_sub"
                    ]
                    for reaction_name in co_bond_reactions_list:
                        if checker.check_reaction(reaction_name, rsmi):
                            co_bond_reaction = True
                            if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            break

                    if co_bond_reaction:
                        # Structural constraint: C-O bond formation at last stage
                        if {"type": "positional", "details": {"target": "C-O bond formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "C-O bond formation", "position": "last_stage"}})

                        # Check if the indole reactant contains a nucleophilic oxygen group
                        nucleophilic_oxygen_fgs = [
                            "Phenol",
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol"
                        ]
                        indole_has_oxygen_nucleophile = False
                        for fg_name in nucleophilic_oxygen_fgs:
                            if checker.check_fg(fg_name, indole_reactant):
                                indole_has_oxygen_nucleophile = True
                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)

                        if indole_has_oxygen_nucleophile:
                            # Check if product has ether group
                            if checker.check_fg("Ether", product_smiles):
                                if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Ether")
                                # Structural constraint: C-O bond formation and Ether co-occurrence
                                if {"type": "co-occurrence", "details": {"targets": ["C-O bond formation", "Ether"]}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["C-O bond formation", "Ether"]}})
                                found_pattern = True
                        else:
                            # Check if indole has a leaving group (electrophile)
                            leaving_group_fgs = [
                                "Primary halide",
                                "Secondary halide",
                                "Tertiary halide",
                                "Aromatic halide",
                                "Triflate",
                                "Mesylate",
                                "Tosylate"
                            ]
                            indole_has_leaving_group = False
                            for fg_name in leaving_group_fgs:
                                if checker.check_fg(fg_name, indole_reactant):
                                    indole_has_leaving_group = True
                                    if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append(fg_name)

                            if indole_has_leaving_group:
                                # Check if another reactant has an oxygen nucleophile
                                for other_reactant in reactants_smiles:
                                    if other_reactant != indole_reactant:
                                        other_reactant_has_oxygen_nucleophile = False
                                        for fg_name in nucleophilic_oxygen_fgs:
                                            if checker.check_fg(fg_name, other_reactant):
                                                other_reactant_has_oxygen_nucleophile = True
                                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                                break

                                        if other_reactant_has_oxygen_nucleophile:
                                            # Check if product has ether group
                                            if checker.check_fg("Ether", product_smiles):
                                                if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                                                    findings_json["atomic_checks"]["functional_groups"].append("Ether")
                                                # Structural constraint: C-O bond formation and Ether co-occurrence
                                                if {"type": "co-occurrence", "details": {"targets": ["C-O bond formation", "Ether"]}} not in findings_json["structural_constraints"]:
                                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["C-O bond formation", "Ether"]}})
                                                found_pattern = True
                                                break

        # Continue traversing
        for i, child in enumerate(node.get("children", [])):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_pattern, findings_json