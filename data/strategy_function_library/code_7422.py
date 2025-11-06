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


# Refactored lists as module-level constants
HETEROCYCLE_SCAFFOLDS = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazine",
    "purine", "quinoline", "isoquinoline", "indole", "benzimidazole",
    "imidazole", "pyrazole", "oxazole", "thiazole", "triazole", "tetrazole",
]

HALOGEN_FGS = ["Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide"]

SUBSTITUTION_REACTIONS_OF_INTEREST = [
    "Aromatic substitution of bromine by chlorine", "Aromatic dehalogenation",
    "N-arylation", "Buchwald-Hartwig", "Suzuki", "Negishi", "Stille",
    "Hiyama-Denmark Coupling", "Sonogashira", "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol", "Ullmann-Goldberg Substitution aryl alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy that starts with a halogenated heterocyclic scaffold
    and elaborates it through sequential substitution reactions.
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

    # Track key features
    has_halogenated_heterocycle_start = False
    substitution_reactions = 0
    maintains_core_scaffold = True
    ring_destruction_occurred = False # New flag for negation constraint

    # Track heterocycle type for each molecule
    molecule_heterocycles = {}

    # Starting materials
    starting_materials = []

    # Track the final product
    final_product = None

    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_halogenated_heterocycle_start, substitution_reactions, maintains_core_scaffold, final_product, findings_json, ring_destruction_occurred

        if path is None:
            path = []

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            if depth == 0:
                final_product = mol_smiles
                print(f"Final product: {mol_smiles}")

            present_heterocycles = []
            for ring in HETEROCYCLE_SCAFFOLDS:
                if checker.check_ring(ring, mol_smiles):
                    present_heterocycles.append(ring)
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)

            if present_heterocycles:
                molecule_heterocycles[mol_smiles] = present_heterocycles
                print(
                    f"Molecule at depth {depth} has heterocycles: {', '.join(present_heterocycles)}"
                )

            if node.get("in_stock", False) or len(node.get("children", [])) == 0:
                has_halogen = False
                found_halogen_fg = None
                for fg in HALOGEN_FGS:
                    if checker.check_fg(fg, mol_smiles):
                        has_halogen = True
                        found_halogen_fg = fg
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                        print(f"Starting material has {fg}: {mol_smiles}")
                        break

                has_heterocycle = len(present_heterocycles) > 0

                if has_halogen and has_heterocycle:
                    print(f"Found halogenated heterocycle starting material: {mol_smiles}")
                    has_halogenated_heterocycle_start = True
                    starting_materials.append(mol_smiles)
                    # Record the co-occurrence constraint if found in a starting material
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "any_halogen_functional_group",
                                "any_heterocyclic_scaffold"
                            ],
                            "scope_description": "Applies to at least one starting material in the route."
                        }
                    })

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                product_heterocycles = []
                for ring in HETEROCYCLE_SCAFFOLDS:
                    if checker.check_ring(ring, product_smiles):
                        product_heterocycles.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                reactant_heterocycles = []
                for r_smiles in reactants_smiles:
                    for ring in HETEROCYCLE_SCAFFOLDS:
                        if checker.check_ring(ring, r_smiles):
                            reactant_heterocycles.append(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                if reactant_heterocycles:
                    print(f"Reactant heterocycles: {', '.join(reactant_heterocycles)}")
                if product_heterocycles:
                    print(f"Product heterocycles: {', '.join(product_heterocycles)}")

                if reactant_heterocycles and not product_heterocycles:
                    print(f"Core scaffold not maintained in reaction at depth {depth}")
                    maintains_core_scaffold = False
                    ring_destruction_occurred = True # Set flag for negation constraint
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                is_substitution = False

                for rxn_type in SUBSTITUTION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction at depth {depth}")
                        is_substitution = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if is_substitution:
                    substitution_reactions += 1
                    print(f"Substitution reaction detected at depth {depth}: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] == "mol": # If current node is chemical, depth increases for reaction child
                new_depth = depth + 1
            # If current node is reaction, depth remains same for chemical child
            
            dfs_traverse(child, new_depth, path + [node])

    dfs_traverse(route)

    strategy_present = (
        has_halogenated_heterocycle_start
        and substitution_reactions >= 1
        and maintains_core_scaffold
    )

    # Add structural constraints based on final flags
    if substitution_reactions >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "substitution_reaction_of_interest",
                "operator": ">=",
                "value": 1
            }
        })
    
    if ring_destruction_occurred:
        # If ring destruction occurred, the negation constraint is NOT met, but we record the finding
        # The overall strategy_present will be False due to maintains_core_scaffold = False
        pass # The 'ring_destruction' atomic check is already added if it occurs
    else:
        # If ring destruction did NOT occur, then the negation constraint is met
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction"
            }
        })

    print(f"Halogenated heterocycle elaboration strategy detected: {strategy_present}")
    print(f"Has halogenated heterocycle start: {has_halogenated_heterocycle_start}")
    print(f"Substitution reactions: {substitution_reactions}")
    print(f"Maintains core scaffold: {maintains_core_scaffold}")

    return strategy_present, findings_json
