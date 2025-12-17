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

HETEROCYCLES_OF_INTEREST = [
    "thiazole",
    "benzothiazole",
    "benzoxazole",
    "benzimidazole",
    "oxazole",
    "isoxazole",
    "imidazole",
    "pyrazole",
    "triazole",
    "tetrazole",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "{thiazole}",
    "{benzothiazole}",
    "{benzothiazole_derivatives_carboxylic-acid/ester}",
    "{benzothiazole_derivatives_aldehyde}",
    "{benzoxazole}",
    "{benzoxazole_arom-aldehyde}",
    "{benzoxazole_carboxylic-acid}",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}",
    "{tetrazole_terminal}",
    "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}",
    "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}",
    "{pyrazole}",
    "{oxadiazole}",
    "{imidazole}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (final 3 steps) formation of specific nitrogen- or sulfur-containing heterocycles. The function identifies formation by either matching a known named reaction (from HETEROCYCLE_FORMATION_REACTIONS) or by structurally confirming the appearance of a heterocycle (from HETEROCYCLES_OF_INTEREST) in the product that was not present in the reactants.
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

    heterocycle_formed = False
    depth_of_formation = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, depth_of_formation, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # First check if this is a known heterocycle formation reaction
                for reaction_name in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        heterocycle_formed = True
                        depth_of_formation = depth
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                # Check for heterocycle appearance by structure
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    # Check if reactants contain the heterocycle
                    reactants_have_heterocycle = False
                    for r_smiles in reactants_smiles.split("."):
                        if checker.check_ring(heterocycle, r_smiles):
                            reactants_have_heterocycle = True
                            break

                    # Check if product contains the heterocycle
                    product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)

                    # If heterocycle is in product but not in reactants, it was formed
                    if product_has_heterocycle and not reactants_have_heterocycle:
                        heterocycle_formed = True
                        depth_of_formation = depth
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                        # Add a generic 'ring_formation' if not already present
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if heterocycle formation occurred in a late stage (within first 3 steps)
    result = heterocycle_formed and depth_of_formation <= 3 and depth_of_formation >= 0

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "heterocycle_formation",
                "position": "late_stage",
                "description": "The formation event must occur at a depth less than or equal to 3 from the final product."
            }
        })

    return result, findings_json
