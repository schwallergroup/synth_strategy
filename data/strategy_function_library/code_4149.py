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
    Detects if the synthesis route uses a convergent approach with C-S bond formation.
    Convergent synthesis means multiple complex fragments are joined in late-stage reactions.
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

    # Track if the synthesis is convergent (multiple complex fragments joined in a late stage)
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json

        # For reaction nodes, check if it's a C-S bond formation
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            is_cs_bond_formation = False
            # Check if this is a C-S bond formation reaction
            if checker.check_reaction("S-alkylation of thiols", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("S-alkylation of thiols")
                is_cs_bond_formation = True
            if checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("S-alkylation of thiols (ethyl)")
                is_cs_bond_formation = True
            if checker.check_reaction("S-alkylation of thiols with alcohols", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("S-alkylation of thiols with alcohols")
                is_cs_bond_formation = True
            if checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("S-alkylation of thiols with alcohols (ethyl)")
                is_cs_bond_formation = True

            product_smiles = rsmi.split(">>", 1)[1]
            reactant_smiles = rsmi.split(">", 1)[0].split(".")

            monosulfide_in_product = checker.check_fg("Monosulfide", product_smiles)
            monosulfide_in_reactants = all(
                checker.check_fg("Monosulfide", r) for r in reactant_smiles
            )

            if monosulfide_in_product:
                if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")

            if monosulfide_in_product and not monosulfide_in_reactants:
                # This implies a C-S bond formation if a monosulfide is formed and not present in all reactants
                if "C-S bond formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("C-S bond formation")
                is_cs_bond_formation = True

            if is_cs_bond_formation:
                # Check if this is a convergent step by examining reactants
                reactants = rsmi.split(">", 1)[0].split(".")
                if len(reactants) >= 2:
                    # If we have at least two reactants and they're both complex (not simple reagents)
                    complex_reactants_count = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol and mol.GetNumHeavyAtoms() > 5
                        ):  # Consider molecules with >5 heavy atoms as complex
                            complex_reactants_count += 1

                    if complex_reactants_count >= 2:
                        # Record the count constraint if met
                        if {"type": "count", "details": {"target": "complex_reactants_in_cs_formation", "operator": ">=", "value": 2, "description": "A C-S bond formation step is considered convergent if it joins at least two complex reactants, defined as molecules with more than 5 heavy atoms."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactants_in_cs_formation", "operator": ">=", "value": 2, "description": "A C-S bond formation step is considered convergent if it joins at least two complex reactants, defined as molecules with more than 5 heavy atoms."}})

                        if depth <= 2:
                            is_convergent = True
                            # Record the positional constraint if met
                            if {"type": "positional", "details": {"target": "convergent_cs_formation", "position": "late_stage", "max_depth": 2, "description": "The convergent C-S bond formation must occur within the last three steps of the synthesis (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "convergent_cs_formation", "position": "late_stage", "max_depth": 2, "description": "The convergent C-S bond formation must occur within the last three steps of the synthesis (depth <= 2)."}})

                            # Record the co-occurrence constraint if all conditions are met
                            if {"type": "co-occurrence", "details": {"targets": ["C-S bond formation", "convergent_step", "late_stage_reaction"], "description": "A single reaction step must be a C-S bond formation, be convergent (joining at least two complex fragments), and occur late in the synthesis (within the last 3 steps)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["C-S bond formation", "convergent_step", "late_stage_reaction"], "description": "A single reaction step must be a C-S bond formation, be convergent (joining at least two complex fragments), and occur late in the synthesis (within the last 3 steps)."}})

        # Recursively traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node.get("type") != "reaction": # If current node is not reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return is_convergent, findings_json
