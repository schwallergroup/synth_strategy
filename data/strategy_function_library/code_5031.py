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
    "indole",
    "quinoline",
    "isoquinoline",
    "purine",
    "thiophene",
    "thiopyran",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a heterocycle from the HETEROCYCLES_OF_INTEREST list is formed in the final step of the synthesis.
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

    final_step_forms_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_forms_heterocycle, findings_json

        # Process current node only if it's the final reaction step (depth=1)
        if node["type"] == "reaction" and depth == 1:
            try:
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]

                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Check for heterocycle formation
                    reactant_smiles = reactants_part.split(".")
                    product_smiles = product_part

                    # Track which heterocycles are present in reactants
                    reactant_heterocycles = set()
                    for r_smiles in reactant_smiles:
                        for het_type in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(het_type, r_smiles):
                                reactant_heterocycles.add(het_type)

                    # Check which heterocycles are in the product
                    product_heterocycles = set()
                    for het_type in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(het_type, product_smiles):
                            product_heterocycles.add(het_type)

                    # Check if any new heterocycles were formed
                    new_heterocycles = product_heterocycles - reactant_heterocycles
                    if new_heterocycles:
                        final_step_forms_heterocycle = True
                        # Record the specific heterocycles formed
                        for het_type in new_heterocycles:
                            if het_type not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(het_type)

            except Exception:
                pass

        # Continue DFS traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for its children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if final_step_forms_heterocycle:
        # Add the structural constraint if a heterocycle was formed in the final step
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "last_stage"
            }
        })

    return final_step_forms_heterocycle, findings_json
