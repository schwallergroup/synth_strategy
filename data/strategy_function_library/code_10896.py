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


# Refactoring for Enumeration: Isolate the list of heterocycles.
HETEROCYCLES_OF_INTEREST = [
    "indole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "imidazole",
    "benzimidazole",
    "quinoline",
    "isoquinoline",
    "thiazole",
    "benzothiazole",
    "oxazole",
    "benzoxazole",
    "triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where two distinct heterocyclic fragments
    are connected via a piperazine linker. The specific heterocycles of interest
    are defined in the HETEROCYCLES_OF_INTEREST list.
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

    # Track if we found the key fragments and connections
    has_piperazine = False
    heterocycle_piperazine_connections = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, heterocycle_piperazine_connections, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperazine
            if checker.check_ring("piperazine", mol_smiles):
                has_piperazine = True
                if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                # Check which heterocycles are in this molecule
                current_heterocycles = set()
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, mol_smiles):
                        current_heterocycles.add(heterocycle)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # If this molecule has piperazine and any heterocycles, track them
                if current_heterocycles:
                    for heterocycle in current_heterocycles:
                        heterocycle_piperazine_connections.add(heterocycle)

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

            # Check for reactions that could connect heterocycles to piperazine
            is_connecting_reaction = False

            # Check for N-arylation reactions
            if checker.check_reaction(
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rxn_smiles
            ):
                is_connecting_reaction = True
                if "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)")

            # Check for alkylation reactions
            if checker.check_reaction("Alkylation of amines", rxn_smiles):
                is_connecting_reaction = True
                if "Alkylation of amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alkylation of amines")

            # Check for reductive amination reactions
            if checker.check_reaction("reductive amination", rxn_smiles):
                is_connecting_reaction = True
                if "reductive amination" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("reductive amination")

            if is_connecting_reaction:
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                has_piperazine_reactant = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )
                if has_piperazine_reactant and "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                has_piperazine_product = checker.check_ring("piperazine", product)
                if has_piperazine_product and "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                heterocycle_in_reactants = set()
                for r in reactants:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, r):
                            heterocycle_in_reactants.add(heterocycle)
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                heterocycle_in_product = set()
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product):
                        heterocycle_in_product.add(heterocycle)
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                # If product has both piperazine and heterocycle, this is a connecting reaction
                if has_piperazine_product and heterocycle_in_product:
                    for heterocycle in heterocycle_in_product:
                        heterocycle_piperazine_connections.add(heterocycle)

                # If a reactant has piperazine and another has a heterocycle, this is a connecting reaction
                elif has_piperazine_reactant and heterocycle_in_reactants:
                    for heterocycle in heterocycle_in_reactants:
                        if has_piperazine_product:  # Only count if the product still has piperazine
                            heterocycle_piperazine_connections.add(heterocycle)

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # The strategy is detected if piperazine is present and connected to at least two different heterocycles.
    strategy_detected = has_piperazine and len(heterocycle_piperazine_connections) >= 2

    if has_piperazine:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "piperazine"
                ]
            }
        })
    
    if len(heterocycle_piperazine_connections) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "distinct_heterocycles_connected_to_piperazine",
                "operator": ">=",
                "value": 2
            }
        })

    return strategy_detected, findings_json
