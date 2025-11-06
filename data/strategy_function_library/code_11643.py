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

HETEROCYCLES_OF_INTEREST = [
    "triazole", "tetrazole", "pyrazole", "imidazole", "oxazole",
    "thiazole", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
    "pyrrole", "furan", "thiophene", "pyridine", "pyrimidine", "pyrazine",
    "pyridazine", "benzimidazole", "benzoxazole", "benzothiazole", "indole",
    "benzofuran", "benzothiophene",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
    "Azide-nitrile click cycloaddition to tetrazole", "Azide-nitrile click cycloaddition to triazole",
    "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}", "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{Huisgen_Ru-catalyzed_1,5_subst}", "{Huisgen_disubst-alkyne}",
    "{1,2,4-triazole_acetohydrazide}", "{1,2,4-triazole_carboxylic-acid/ester}",
    "{pyrazole}", "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}", "{benzothiazole}",
    "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{thiazole}", "{imidazole}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a convergent synthesis where a heterocycle is formed
    in the final step by combining two separately synthesized fragments.
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

    found_convergent_synthesis = False
    found_heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_synthesis, found_heterocycle_formation, findings_json

        # Check the final reaction step
        if node["type"] == "reaction" and depth == 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if product contains a heterocycle that wasn't in the reactants
                product_heterocycles = []
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.append(heterocycle)

                reactant_heterocycles = set()
                for reactant in reactants:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)

                # Check if there's a heterocycle in the product that wasn't in any reactant
                new_heterocycles = [
                    h for h in product_heterocycles if h not in reactant_heterocycles
                ]
                if new_heterocycles:
                    found_heterocycle_formation = True
                    # Record atomic checks for newly formed heterocycles
                    for h in new_heterocycles:
                        if h not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(h)
                    # Record positional constraint for ring_formation at last_stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "position": "last_stage"
                        }
                    })

                # Check if this is a convergent synthesis (multiple significant branches)
                if len(node.get("children", [])) >= 2:
                    significant_branches = 0
                    for child in node.get("children", []):
                        # A significant branch has at least one reaction step
                        if child.get("type") == "mol" and len(child.get("children", [])) > 0:
                            significant_branches += 1

                    if significant_branches >= 2:
                        found_convergent_synthesis = True
                        # Record count constraint for significant_reactant_branches
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "significant_reactant_branches",
                                "operator": ">=",
                                "value": 2
                            }
                        })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'chemical' (mol), depth increases
                new_depth = depth + 1
            # If current node is 'reaction', depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    result = found_convergent_synthesis and found_heterocycle_formation

    # Record co-occurrence constraint if both conditions are met
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "final_step_heterocycle_formation",
                    "convergent_final_step"
                ]
            }
        })

    return result, findings_json
