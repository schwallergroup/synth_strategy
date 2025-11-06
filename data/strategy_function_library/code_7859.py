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
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane",
    "dioxathiolane", "thiazolidine", "oxazolidine", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole", "benzoxazole",
    "benzothiazole", "benzimidazole", "pteridin", "phenothiazine",
    "phenoxazine", "dibenzofuran", "dibenzothiophene", "xanthene",
    "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine",
    "porphyrin", "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of a new heterocycle in a single reaction step during the early stages of a synthesis (depth >= 3). The check is performed against a predefined list of heterocycle structures specified in HETEROCYCLES_OF_INTEREST.
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

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, findings_json

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains heterocycles not present in reactants
                product_heterocycles = set()
                reactants_heterocycles = set()

                # Check heterocycles in product
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.add(heterocycle)
                        # Add to findings if it's a new heterocycle found in product
                        # This will be filtered later by new_heterocycles

                # Check heterocycles in reactants
                for reactant in reactants_smiles:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, reactant):
                            reactants_heterocycles.add(heterocycle)

                # Check if new heterocycles were formed
                new_heterocycles = product_heterocycles - reactants_heterocycles
                if new_heterocycles:
                    heterocycle_formed = True
                    for h in new_heterocycles:
                        if h not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(h)
                    # Add structural constraint if a new heterocycle is formed at early stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "position": "early_stage (depth >= 3)"
                        }
                    })

        for child in node.get("children", []):
            if heterocycle_formed: # Optimization to stop searching once found
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return heterocycle_formed, findings_json
