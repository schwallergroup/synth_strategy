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

HETEROCYCLES_OF_INTEREST = [
    "furan",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "indole",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
    "thiophene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy to connect heterocycles via a thioether bridge, formed by an S-alkylation reaction. The heterocycles of interest include: furan, pyrrole, pyridine, pyrazole, imidazole, oxazole, thiazole, pyrimidine, pyrazine, pyridazine, triazole, tetrazole, indole, benzimidazole, benzoxazole, benzothiazole, quinoline, isoquinoline, and thiophene.
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

    thioether_connection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal thioether_connection_found, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if checker.check_fg("Monosulfide", product):
                if "Monosulfide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Monosulfide")

                heterocycle_count = 0
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product):
                        heterocycle_count += 1
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                if heterocycle_count >= 2:
                    # Structural Constraint: heterocycle_of_interest_in_product >= 2
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "heterocycle_of_interest_in_product",
                            "operator": ">=",
                            "value": 2
                        }
                    })

                is_thioether_formation = False
                s_alkylation_reactions = [
                    "S-alkylation of thiols",
                    "S-alkylation of thiols (ethyl)",
                    "S-alkylation of thiols with alcohols",
                    "S-alkylation of thiols with alcohols (ethyl)"
                ]
                for rxn_name in s_alkylation_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_thioether_formation = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)

                has_thioether_in_reactants = False
                for r in reactants:
                    if checker.check_fg("Monosulfide", r):
                        has_thioether_in_reactants = True
                        break
                
                if has_thioether_in_reactants == False:
                    # Structural Constraint: negation of Monosulfide in reactants
                    findings_json["structural_constraints"].append({
                        "type": "negation",
                        "details": {
                            "target": "Monosulfide",
                            "scope": "reactants"
                        }
                    })

                reactant_heterocycles = []
                for r in reactants:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, r):
                            reactant_heterocycles.append((r, heterocycle))
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                if len(reactant_heterocycles) >= 1:
                    # Structural Constraint: heterocycle_of_interest_in_reactant >= 1
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "heterocycle_of_interest_in_reactant",
                            "operator": ">=",
                            "value": 1
                        }
                    })

                if (
                    is_thioether_formation
                    and len(reactant_heterocycles) >= 1
                    and not has_thioether_in_reactants
                ):
                    thioether_connection_found = True

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # This means it's a chemical node
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return thioether_connection_found, findings_json
