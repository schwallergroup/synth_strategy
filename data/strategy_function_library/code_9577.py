from typing import Tuple, Dict, List
import copy
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


HETEROCYCLES_OF_INTEREST = [
    "benzofuran", "imidazole", "benzimidazole", "oxazole", "thiazole",
    "furan", "pyrrole", "pyridine", "pyrazole", "thiophene", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole", "triazole", "tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving the alkylation of a polyol scaffold (a molecule with two or more alcohol groups) with a specific heterocyclic fragment, followed by a later-stage carboxylic acid protection. The strategy specifically looks for heterocycles from the following list: benzofuran, imidazole, benzimidazole, oxazole, thiazole, furan, pyrrole, pyridine, pyrazole, thiophene, isoxazole, isothiazole, oxadiazole, thiadiazole, triazole, tetrazole.
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

    # Track if we found the key components
    found_diol = False
    found_heterocycle = False
    found_heterocycle_name = None
    found_alkylation = False
    found_ester_protection = False

    # Track depths for sequence verification
    alkylation_depth = float("inf")
    protection_depth = float("inf")

    # Track molecules for verification
    diol_smiles = None
    heterocycle_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal found_diol, found_heterocycle, found_alkylation, found_ester_protection
        nonlocal alkylation_depth, protection_depth, diol_smiles, heterocycle_smiles, found_heterocycle_name
        nonlocal findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for polyol scaffold (molecule with multiple alcohol groups)
            if not found_diol:
                # Look for molecules with multiple alcohol groups
                alcohol_groups = ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol"]
                alcohol_count = 0
                detected_alcohols = []
                for fg in alcohol_groups:
                    if checker.check_fg(fg, mol_smiles):
                        alcohol_count += 1
                        detected_alcohols.append(fg)

                if alcohol_count >= 2:
                    found_diol = True
                    diol_smiles = mol_smiles
                    findings_json["atomic_checks"]["functional_groups"].extend(detected_alcohols)
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": ["Primary alcohol", "Secondary alcohol", "Tertiary alcohol"], "operator": ">=", "value": 2, "scope": "molecule", "description": "A molecule with two or more alcohol groups (polyol) must be present in the route."}})
                    print(f"Found polyol scaffold at depth {depth}: {mol_smiles}")

            # Check for heterocyclic fragments
            if not found_heterocycle:
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, mol_smiles):
                        found_heterocycle = True
                        found_heterocycle_name = ring
                        heterocycle_smiles = mol_smiles
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                        print(f"Found heterocycle ({ring}) at depth {depth}: {mol_smiles}")
                        break

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkylation reaction (Williamson ether synthesis or similar)
                if not found_alkylation:
                    # Check for Williamson ether synthesis
                    if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                        found_alkylation = True
                        alkylation_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                        print(f"Found Williamson ether synthesis at depth {depth}")
                    # Check for other ether formation reactions
                    elif any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    ) and checker.check_fg("Ether", product):
                        found_alkylation = True
                        alkylation_depth = depth
                        for r in reactants:
                            if checker.check_fg("Primary alcohol", r):
                                findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                            if checker.check_fg("Secondary alcohol", r):
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                            if checker.check_fg("Tertiary alcohol", r):
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")
                        print(f"Found alkylation reaction (C-O bond formation) at depth {depth}")

                # Check for ester protection or esterification
                if not found_ester_protection:
                    # Check for specific protection/esterification reactions
                    if (
                        checker.check_reaction("Protection of carboxylic acid", rsmi)
                    ):
                        found_ester_protection = True
                        protection_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("Protection of carboxylic acid")
                        print(f"Found carboxylic acid protection/esterification at depth {depth}")
                    elif (
                        checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    ):
                        found_ester_protection = True
                        protection_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("Esterification of Carboxylic Acids")
                        print(f"Found carboxylic acid protection/esterification at depth {depth}")
                    elif (
                        checker.check_reaction("O-alkylation of carboxylic acids with diazo compounds", rsmi)
                    ):
                        found_ester_protection = True
                        protection_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("O-alkylation of carboxylic acids with diazo compounds")
                        print(f"Found carboxylic acid protection/esterification at depth {depth}")
                    # Check for general ester formation
                    elif any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    ) and checker.check_fg("Ester", product):
                        found_ester_protection = True
                        protection_depth = depth
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        print(f"Found ester formation at depth {depth}")
                    # Check if the product already contains an ester group
                    elif checker.check_fg("Ester", product) and not all(
                        checker.check_fg("Ester", r) for r in reactants
                    ):
                        found_ester_protection = True
                        protection_depth = depth
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")
                        print(f"Found ester group in product at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all components are found
    all_components_found = found_diol and found_heterocycle and found_alkylation

    # "followed by" means protection happens at a later stage, which is a SMALLER depth value.
    # e.g., alkylation at depth 5, protection at depth 3. So protection_depth < alkylation_depth
    correct_sequence = True
    if found_alkylation and found_ester_protection:
        correct_sequence = protection_depth < alkylation_depth
        if correct_sequence:
            findings_json["structural_constraints"].append({"type": "sequence", "details": {"before": "alkylation_reaction", "after": "ester_protection_reaction", "description": "An ester protection/formation reaction must occur at a later synthetic stage (smaller depth) than the alkylation reaction."}})

    if all_components_found:
        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["polyol_scaffold", "heterocycle_of_interest", "alkylation_reaction"], "scope": "route", "description": "The route must contain a polyol scaffold, one of the specified heterocycles, and an alkylation reaction."}})

    print(f"Found diol: {found_diol}")
    print(f"Found heterocycle: {found_heterocycle} ({found_heterocycle_name})")
    print(f"Found alkylation: {found_alkylation} (depth: {alkylation_depth})")
    print(f"Found ester protection: {found_ester_protection} (depth: {protection_depth})")
    print(f"Correct sequence: {correct_sequence}")

    # Return True if all key components of the strategy are found
    result = all_components_found and correct_sequence
    return result, findings_json
