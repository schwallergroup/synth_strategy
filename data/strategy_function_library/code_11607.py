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


ALKYLATING_FGS = [
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
    "Triflate",
    "Mesylate",
    "Tosylate",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the formation of a quaternary carbon center alpha to a nitrile group.
    It identifies reactions where a nitrile-containing reactant is treated with a common alkylating agent,
    resulting in the formation of a new quaternary center. The function checks for the presence of specific
    alkylating functional groups, including: Primary halide, Secondary halide, Tertiary halide, Triflate,
    Mesylate, and Tosylate.
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

    quaternary_center_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal quaternary_center_formed, findings_json

        if node["type"] == "reaction" and not quaternary_center_formed:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile in both reactants and product
                reactant_has_nitrile = False
                for r in reactants:
                    if checker.check_fg("Nitrile", r):
                        reactant_has_nitrile = True
                        if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                product_has_nitrile = checker.check_fg("Nitrile", product)
                if product_has_nitrile:
                    if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

                # Check for potential alkylating agents in reactants
                alkylating_agent_present = False
                for fg in ALKYLATING_FGS:
                    for r in reactants:
                        if checker.check_fg(fg, r):
                            alkylating_agent_present = True
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                            break
                    if alkylating_agent_present: # Optimization: if found, no need to check other reactants for this FG
                        break

                # If we have nitrile in both reactants and product and potential alkylation conditions
                if (
                    reactant_has_nitrile
                    and product_has_nitrile
                    and alkylating_agent_present
                ):
                    try:
                        # Convert to RDKit molecules for detailed analysis
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [
                            Chem.MolFromSmiles(r)
                            for r in reactants
                            if checker.check_fg("Nitrile", r)
                        ]

                        if not reactant_mols:
                            return

                        # Find nitrile carbon atoms in the product
                        product_nitrile_carbons = []
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "C":
                                for neighbor in atom.GetNeighbors():
                                    if neighbor.GetSymbol() == "N" and neighbor.GetDegree() == 1:
                                        # This carbon is part of a nitrile group
                                        product_nitrile_carbons.append(atom.GetIdx())

                        # For each nitrile carbon in the product
                        for nitrile_carbon_idx in product_nitrile_carbons:
                            nitrile_carbon = product_mol.GetAtomWithIdx(nitrile_carbon_idx)

                            # Find alpha carbons (directly attached to nitrile carbon)
                            alpha_carbons = []
                            for neighbor in nitrile_carbon.GetNeighbors():
                                if neighbor.GetSymbol() == "C":
                                    alpha_carbons.append(neighbor)

                            # Check if any alpha carbon is quaternary (has 4 neighbors)
                            for alpha_carbon in alpha_carbons:
                                if alpha_carbon.GetDegree() == 4:
                                    # Now check if this alpha carbon was tertiary in any reactant
                                    alpha_carbon_quaternary_in_reactant = False

                                    # For each reactant with a nitrile
                                    for reactant_mol in reactant_mols:
                                        reactant_nitrile_carbons = []
                                        for atom in reactant_mol.GetAtoms():
                                            if atom.GetSymbol() == "C":
                                                for neighbor in atom.GetNeighbors():
                                                    if (
                                                        neighbor.GetSymbol() == "N"
                                                        and neighbor.GetDegree() == 1
                                                    ):
                                                        reactant_nitrile_carbons.append(
                                                            atom.GetIdx()
                                                        )

                                        # For each nitrile carbon in this reactant
                                        for r_nitrile_carbon_idx in reactant_nitrile_carbons:
                                            r_nitrile_carbon = reactant_mol.GetAtomWithIdx(
                                                r_nitrile_carbon_idx
                                            )

                                            # Find alpha carbons in reactant
                                            r_alpha_carbons = []
                                            for neighbor in r_nitrile_carbon.GetNeighbors():
                                                if neighbor.GetSymbol() == "C":
                                                    r_alpha_carbons.append(neighbor)

                                            # Check if any alpha carbon is already quaternary
                                            for r_alpha_carbon in r_alpha_carbons:
                                                if r_alpha_carbon.GetDegree() == 4:
                                                    alpha_carbon_quaternary_in_reactant = True
                                                    break

                                    # If alpha carbon was not quaternary in reactants but is in product
                                    if not alpha_carbon_quaternary_in_reactant:
                                        print(
                                            f"Quaternary center formation via nitrile alpha-carbon alkylation detected"
                                        )
                                        quaternary_center_formed = True
                                        if "quaternary_alpha_nitrile_alkylation" not in findings_json["atomic_checks"]["named_reactions"]:
                                            findings_json["atomic_checks"]["named_reactions"].append("quaternary_alpha_nitrile_alkylation")
                                        # Add structural constraint if detected
                                        findings_json["structural_constraints"].append({
                                            "type": "count",
                                            "details": {
                                                "target": "quaternary_alpha_nitrile_alkylation",
                                                "operator": ">=",
                                                "value": 1
                                            }
                                        })
                                        return
                    except Exception as e:
                        print(f"Error analyzing molecule structure: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for children
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return quaternary_center_formed, findings_json
