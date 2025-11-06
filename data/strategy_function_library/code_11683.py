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


HETEROCYCLE_TYPES_FOR_CONVERGENCE = [
    "furan",
    "pyran",
    "dioxane",
    "tetrahydrofuran",
    "tetrahydropyran",
    "oxirane",
    "oxetane",
    "oxolane",
    "oxane",
    "dioxolane",
    "dioxolene",
    "trioxane",
    "dioxepane",
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
    "pyrrolidine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "aziridine",
    "azetidine",
    "azepane",
    "diazepane",
    "indole",
    "quinoline",
    "isoquinoline",
    "purine",
    "carbazole",
    "acridine",
    "thiophene",
    "thiopyran",
    "thiirane",
    "thietane",
    "thiolane",
    "thiane",
    "dithiane",
    "dithiolane",
    "benzothiophene",
    "oxathiolane",
    "dioxathiolane",
    "thiazolidine",
    "oxazolidine",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "indazole",
    "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a convergent synthesis strategy involving multiple heterocyclic fragments.
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

    # Track heterocycles found in each reaction and overall
    heterocycle_fragments = set()
    convergent_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_fragments, convergent_reactions, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Track heterocycles in this reaction
                reaction_heterocycles = set()

                # Check reactants for heterocycles
                for reactant in reactants:
                    for heterocycle in HETEROCYCLE_TYPES_FOR_CONVERGENCE:
                        if checker.check_ring(heterocycle, reactant):
                            reaction_heterocycles.add(heterocycle)
                            heterocycle_fragments.add(heterocycle)
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # Check if this reaction combines multiple heterocycles
                if len(reaction_heterocycles) >= 2:
                    # Verify product contains at least one of the heterocycles
                    product_has_heterocycle = False
                    for heterocycle in reaction_heterocycles:
                        if checker.check_ring(heterocycle, product):
                            product_has_heterocycle = True
                            break

                    if product_has_heterocycle:
                        convergent_reactions += 1

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

    # Count unique heterocyclic fragments
    fragment_count = len(heterocycle_fragments)
    print(f"Found {fragment_count} different heterocyclic fragments: {heterocycle_fragments}")
    print(f"Found {convergent_reactions} convergent reactions")

    # Determine the final result
    result = fragment_count >= 3 and convergent_reactions >= 1

    # Add structural constraints to findings_json if met
    if fragment_count >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_heterocyclic_fragments",
                "operator": ">=",
                "value": 3
            }
        })
    if convergent_reactions >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_reactions",
                "operator": ">=",
                "value": 1
            }
        })

    # Return True if we found evidence of convergent synthesis:
    # At least 3 different heterocycles AND at least 1 convergent reaction
    return result, findings_json
