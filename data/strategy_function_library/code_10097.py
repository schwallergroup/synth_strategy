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


HETEROCYCLIC_RINGS = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "benzoxazole", "benzothiazole", "benzimidazole", "pteridin",
    "phenothiazine", "phenoxazine", "dibenzofuran", "dibenzothiophene",
    "xanthene", "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine",
    "porphyrin", "indazole", "benzotriazole",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "Niementowski_quinazoline", "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "Huisgen_Cu-catalyzed_1,4-subst", "Huisgen_Ru-catalyzed_1,5_subst",
    "1,2,4-triazole_acetohydrazide", "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine", "pyrazole", "Paal-Knorr pyrrole", "triaryl-imidazole",
    "Fischer indole", "Friedlaender chinoline", "benzofuran", "benzothiophene",
    "indole", "oxadiazole", "imidazole", "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition", "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation", "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear synthetic strategy that culminates in the formation of a new heterocycle in the final step. This check is specific to the `HETEROCYCLIC_RINGS` list.
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

    # Track synthesis characteristics
    is_linear = True  # Assume linear until proven otherwise
    has_final_heterocyclization = False
    reaction_count = 0
    final_step_found = False

    def is_final_step(node):
        """Determine if this is the final synthetic step (depth 0)"""
        if "depth" in node["metadata"]:
            return node["metadata"]["depth"] == 0

        # Alternative method: check if the product doesn't appear as a reactant in any other step
        product_smiles = node["metadata"]["mapped_reaction_smiles"].split(">")[-1]

        # Traverse the route to see if this product is used as a reactant elsewhere
        def check_if_used_as_reactant(check_node, target_smiles):
            if check_node["type"] == "reaction" and check_node != node:
                reactants_part = check_node["metadata"]["mapped_reaction_smiles"].split(">")[0]
                if target_smiles in reactants_part.split("."):
                    return True

            for child in check_node.get("children", []):
                if check_if_used_as_reactant(child, target_smiles):
                    return True

            return False

        # If the product is not used as a reactant elsewhere, it's the final step
        return not check_if_used_as_reactant(route, product_smiles)

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, has_final_heterocyclization, reaction_count, final_step_found, findings_json

        if node["type"] == "reaction":
            reaction_count += 1
            print(f"Analyzing reaction at depth {depth}")

            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            is_final = is_final_step(node)
            if is_final:
                print(f"Found final step: {rsmi}")
                final_step_found = True
                # Record positional constraint: ring_formation at last_stage
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "ring_formation",
                        "position": "last_stage"
                    }
                })

            reactant_smiles_list = reactants_part.split(".")
            significant_reactants = []

            for r_smiles in reactant_smiles_list:
                if not r_smiles:
                    continue

                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol is None:
                    print(f"Warning: Could not parse reactant SMILES: {r_smiles}")
                    continue

                if r_mol.GetNumAtoms() > 5 and any(
                    atom.GetSymbol() == "C" for atom in r_mol.GetAtoms()
                ):
                    significant_reactants.append(r_smiles)

            if len(significant_reactants) > 1:
                print(
                    f"Found convergent step with {len(significant_reactants)} significant reactants"
                )
                is_linear = False
                # Record negation constraint: convergent_step
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "convergent_step"
                    }
                })

            if is_final:
                product_mol = Chem.MolFromSmiles(product_part)
                if product_mol is None:
                    print(f"Warning: Could not parse product SMILES: {product_part}")
                else:
                    product_heterocycles = []
                    for ring_name in HETEROCYCLIC_RINGS:
                        if checker.check_ring(ring_name, product_part):
                            product_heterocycles.append(ring_name)
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)

                    if product_heterocycles:
                        print(f"Product contains heterocycles: {', '.join(product_heterocycles)}")

                    reactant_heterocycles = set()
                    for r_smiles in reactant_smiles_list:
                        if not r_smiles:
                            continue

                        for ring_name in product_heterocycles:
                            if checker.check_ring(ring_name, r_smiles):
                                reactant_heterocycles.add(ring_name)

                    new_heterocycles = set(product_heterocycles) - reactant_heterocycles
                    if new_heterocycles:
                        print(f"New heterocycles formed: {', '.join(new_heterocycles)}")
                        has_final_heterocyclization = True
                        # Record named reaction: ring_formation
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (which are chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (which are reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = is_linear and has_final_heterocyclization and reaction_count > 1 and final_step_found
    print(f"Linear synthesis: {is_linear}")
    print(f"Has final heterocyclization: {has_final_heterocyclization}")
    print(f"Multiple reactions: {reaction_count > 1}")
    print(f"Final step found: {final_step_found}")
    print(f"Overall result: {result}")

    # Record count constraint: reaction > 1
    if reaction_count > 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">",
                "value": 1
            }
        })

    return result, findings_json
