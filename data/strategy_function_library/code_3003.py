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
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "indole", "quinoline", "isoquinoline", "thiophene", "benzoxazole",
    "benzothiazole", "benzimidazole",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
    "Benzothiazole formation from aldehyde", "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid", "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide", "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)", "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide", "Benzimidazole formation from ester/carboxylic acid",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the late-stage formation of specific heterocycles. This is determined by
    checking if the final reaction step (depth=1) is a known heterocycle-forming reaction
    type (from HETEROCYCLE_FORMING_REACTIONS) or if it results in the formation of a new
    heterocycle (from HETEROCYCLES_OF_INTEREST) that was not present in the reactants.
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

    has_late_stage_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_heterocycle, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:

            # Check if this is the final reaction (depth 1 in retrosynthetic direction)
            if depth == 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                # First check if this is a known heterocycle-forming reaction
                for rxn_type in HETEROCYCLE_FORMING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_late_stage_heterocycle = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        if {"type": "positional", "details": {"target": "HETEROCYCLE_FORMING_REACTIONS", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "HETEROCYCLE_FORMING_REACTIONS", "position": "last_stage"}})
                        return

                # Parse reactants and product
                reactants = reactants_part.split(".")
                product = products_part

                # Check if any heterocycle is present in the product but not in any reactant
                product_heterocycles = []
                reactant_heterocycles = []

                # Check heterocycles in product
                for cycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(cycle, product):
                        product_heterocycles.append(cycle)
                        if cycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(cycle)

                # Check heterocycles in reactants
                for reactant in reactants:
                    for cycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(cycle, reactant):
                            reactant_heterocycles.append(cycle)
                            # Only add to findings if it's a new finding, not if it's just present in reactant
                            # The goal is to find *newly formed* heterocycles

                # Find heterocycles in product that aren't in any reactant
                new_heterocycles = [
                    cycle for cycle in product_heterocycles if cycle not in reactant_heterocycles
                ]

                if new_heterocycles:
                    has_late_stage_heterocycle = True
                    if {"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}})
                    # The individual new heterocycles are already added to ring_systems in product_heterocycles loop
                    return

        # Continue traversing with modified depth logic
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_heterocycle, findings_json