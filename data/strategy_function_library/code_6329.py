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


HETEROCYCLE_TYPES = [
    "pyrazole",
    "pyrimidine",
    "pyridine",
    "furan",
    "thiophene",
    "imidazole",
    "oxazole",
    "thiazole",
    "triazole",
    "tetrazole",
    "morpholine",
    "piperidine",
    "piperazine",
    "indole",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "pyrrole",
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
    "pyrrolidine",
    "aziridine",
    "azetidine",
    "azepane",
    "diazepane",
    "quinoline",
    "isoquinoline",
    "purine",
]

COUPLING_REACTIONS = [
    "Suzuki",
    "Buchwald-Hartwig",
    "N-arylation",
    "Sonogashira",
    "Heck",
    "Negishi",
    "Stille",
    "Ullmann",
    "Acylation of Nitrogen Nucleophiles",
    "Schotten-Baumann",
    "Mitsunobu",
    "Williamson Ether Synthesis",
    "Amidation",
    "Esterification",
    "Urea synthesis",
    "Thiourea synthesis",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent strategy that couples pre-formed heterocyclic fragments. The function identifies reactions from a specific list of coupling reactions (e.g., Suzuki, Amidation) that involve reactants containing heterocycles from a specific list (e.g., pyridine, furan, indole). A synthesis is flagged if these coupling events outnumber de novo heterocycle formations.
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

    # Track heterocycle construction vs coupling
    heterocycle_construction_count = 0
    fragment_coupling_count = 0
    heterocycles_in_fragments = set()

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_construction_count, fragment_coupling_count, heterocycles_in_fragments, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Identify heterocycles in reactants
                heterocycles_in_reactants = set()
                for reactant_smiles in reactants_smiles:
                    for heterocycle in HETEROCYCLE_TYPES:
                        if checker.check_ring(heterocycle, reactant_smiles):
                            heterocycles_in_reactants.add(heterocycle)
                            heterocycles_in_fragments.add(heterocycle)
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # Identify heterocycles in product
                heterocycles_in_product = set()
                for heterocycle in HETEROCYCLE_TYPES:
                    if checker.check_ring(heterocycle, product_smiles):
                        heterocycles_in_product.add(heterocycle)
                        # If a heterocycle is formed de novo, record it as a ring_formation
                        if heterocycle not in heterocycles_in_reactants and "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check if this is a heterocycle construction
                new_heterocycles = heterocycles_in_product - heterocycles_in_reactants
                if new_heterocycles:
                    heterocycle_construction_count += 1

                # Check for coupling reactions between fragments
                if heterocycles_in_reactants and len(reactants_smiles) > 1:
                    # Check if this is a known coupling reaction
                    is_coupling = False
                    for rxn in COUPLING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_coupling = True
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                    if is_coupling:
                        fragment_coupling_count += 1

            except Exception as e:
                print(f"Error processing reaction: {e}")

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

    # Strategy is detected if:
    # 1. We have more fragment couplings than heterocycle constructions
    # 2. At least 2 different heterocycles are present in fragments
    result = (
        fragment_coupling_count > heterocycle_construction_count
        and len(heterocycles_in_fragments) >= 2
    )

    if fragment_coupling_count > heterocycle_construction_count:
        findings_json["structural_constraints"].append({
            "type": "relative_count",
            "details": {
                "target_a": "fragment_coupling_with_heterocycle",
                "operator": ">",
                "target_b": "de_novo_heterocycle_formation"
            }
        })
    if len(heterocycles_in_fragments) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "unique_reactant_heterocycles",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
