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
    "pteridin",
    "phenothiazine",
    "phenoxazine",
    "dibenzofuran",
    "dibenzothiophene",
    "xanthene",
    "thioxanthene",
    "pyrroline",
    "pyrrolidone",
    "imidazolidine",
    "porphyrin",
    "indazole",
    "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the late-stage formation of specific heterocycles defined in the HETEROCYCLES_OF_INTEREST list. A formation event is flagged if a heterocycle from the list is present in the product but not in any of the reactants. The strategy is considered successful if this occurs within the first two steps of the synthesis (depth <= 2).
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

    heterocycle_connected = False
    depth_of_connection = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_connected, depth_of_connection, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains any heterocycle
                product_mol = Chem.MolFromSmiles(product)
                if product_mol is None:
                    print(f"Warning: Could not parse product SMILES: {product}")
                    return

                # Check for heterocycles in product
                product_heterocycles = []
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.append(heterocycle)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                if product_heterocycles:
                    # Check if any of these heterocycles are not present in all reactants
                    for heterocycle in product_heterocycles:
                        heterocycle_in_any_reactant = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol is None:
                                print(f"Warning: Could not parse reactant SMILES: {reactant}")
                                continue

                            if checker.check_ring(heterocycle, reactant):
                                heterocycle_in_any_reactant = True
                                break

                        # If heterocycle is in product but not in any reactant, it was formed in this step
                        if not heterocycle_in_any_reactant:
                            print(f"Heterocycle {heterocycle} connection detected at depth {depth}")
                            heterocycle_connected = True
                            depth_of_connection = min(depth_of_connection, depth)
                            # Record the reaction as a 'ring_formation' if it leads to a new heterocycle
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Consider it late stage if it happens in the first half of the synthesis
    # (remember lower depth means later in synthesis)
    result = heterocycle_connected and depth_of_connection <= 2
    print(f"Heterocycle connected: {heterocycle_connected}, Depth: {depth_of_connection}")

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "depth <= 2"
            }
        })

    return result, findings_json
