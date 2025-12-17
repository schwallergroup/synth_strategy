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


N_HETEROCYCLES_FOR_ALKYLATION = [
    "pyrazole", "imidazole", "triazole", "tetrazole", "pyrrole", "indole",
    "pyridine", "piperidine", "piperazine", "morpholine", "pyrrolidine",
    "azetidine", "aziridine", "azepane", "diazepane", "benzimidazole",
    "benzotriazole", "indazole", "quinoline", "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects N-alkylation of specific heterocycles using a mesylate as an activating group for an alcohol.
    The strategy is flagged if a mesylate is formed from an alcohol and subsequently used for N-alkylation,
    or if a mesylate is present as a reactant in an N-alkylation step. The specific heterocycles checked
    are defined in the N_HETEROCYCLES_FOR_ALKYLATION list.
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

    # Initialize flags to track strategy components
    mesylate_formation_found = False
    n_alkylation_found = False

    # Track mesylate molecules
    mesylate_molecules = set()

    # Helper function to check if a molecule contains a mesylate group
    def is_mesylate(mol_smiles):
        if checker.check_fg("Mesylate", mol_smiles):
            if "Mesylate" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Mesylate")
            return True
        return False

    # Helper function to check if a molecule contains a heterocycle with N
    def has_n_heterocycle(mol_smiles):
        for ring in N_HETEROCYCLES_FOR_ALKYLATION:
            if checker.check_ring(ring, mol_smiles):
                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                return True
        return False

    # DFS traversal to analyze the route
    def dfs_traverse(node, depth=0):
        nonlocal mesylate_formation_found, n_alkylation_found
        nonlocal mesylate_molecules, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for mesylate formation (alcohol to mesylate)
            alcohol_found_in_reactants = False
            for r in reactants:
                if checker.check_fg("Alcohol", r):
                    alcohol_found_in_reactants = True
                    if "Alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Alcohol")
                    break

            if alcohol_found_in_reactants and is_mesylate(product):
                sulfonic_ester_formation = False
                if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                    sulfonic_ester_formation = True
                    if "Formation of Sulfonic Esters" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Formation of Sulfonic Esters")
                
                tms_sulfonic_ester_formation = False
                if checker.check_reaction("Formation of Sulfonic Esters on TMS protected alcohol", rsmi):
                    tms_sulfonic_ester_formation = True
                    if "Formation of Sulfonic Esters on TMS protected alcohol" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Formation of Sulfonic Esters on TMS protected alcohol")

                if sulfonic_ester_formation or tms_sulfonic_ester_formation:
                    mesylate_formation_found = True
                    canonical_product = Chem.CanonSmiles(product)
                    mesylate_molecules.add(canonical_product)

            # Check for N-alkylation with mesylate
            mesylate_in_reactants = any(is_mesylate(r) for r in reactants)
            n_heterocycle_in_product = has_n_heterocycle(product)

            if mesylate_in_reactants and n_heterocycle_in_product:
                alkylation_of_amines = False
                if checker.check_reaction("Alkylation of amines", rsmi):
                    alkylation_of_amines = True
                    if "Alkylation of amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Alkylation of amines")

                n_heterocycle_in_reactants = any(has_n_heterocycle(r) for r in reactants)

                if (alkylation_of_amines or (n_heterocycle_in_reactants and not all(has_n_heterocycle(r) for r in reactants))):
                    n_alkylation_found = True
                    # Add structural constraint if this condition is met
                    # This constraint is for 'Alkylation of amines' and 'Mesylate' co-occurrence
                    constraint_obj = {
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "Alkylation of amines",
                                "Mesylate"
                            ]
                        }
                    }
                    if constraint_obj not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_obj)

        elif node["type"] == "mol":
            mol_smiles = node["smiles"]
            if is_mesylate(mol_smiles):
                canonical_smiles = Chem.CanonSmiles(mol_smiles)
                mesylate_molecules.add(canonical_smiles)

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (i.e., 'chemical'), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    has_mesylates = len(mesylate_molecules) > 0

    combined_strategy_detected = (mesylate_formation_found and n_alkylation_found) or (
        has_mesylates and n_alkylation_found
    )

    return combined_strategy_detected, findings_json
