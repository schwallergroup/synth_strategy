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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthesis strategy that builds a pyridine scaffold
    and connects it to a piperidine via C-N bond formation.
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
    result = False

    # Track key components and reactions
    pyridine_nodes = []
    piperidine_nodes = []
    connection_nodes = []
    connection_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for pyridine
            if checker.check_ring("pyridine", mol_smiles):
                pyridine_nodes.append((mol_smiles, depth))
                findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                print(f"Found pyridine at depth {depth}: {mol_smiles}")

            # Check for piperidine
            if checker.check_ring("piperidine", mol_smiles):
                piperidine_nodes.append((mol_smiles, depth))
                findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                print(f"Found piperidine at depth {depth}: {mol_smiles}")

            # Check for connected structure (both rings in same molecule)
            if checker.check_ring("pyridine", mol_smiles) and checker.check_ring(
                "piperidine", mol_smiles
            ):
                connection_nodes.append((mol_smiles, depth))
                print(f"Found connected structure at depth {depth}: {mol_smiles}")

        # Process reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-N bond forming reactions
                if checker.check_reaction("N-arylation", rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append("N-arylation")

                    # Check if this reaction connects pyridine and piperidine
                    has_pyridine_reactant = any(
                        checker.check_ring("pyridine", r) for r in reactants
                    )
                    has_piperidine_reactant = any(
                        checker.check_ring("piperidine", r) for r in reactants
                    )
                    has_both_in_product = checker.check_ring(
                        "pyridine", product
                    ) and checker.check_ring("piperidine", product)

                    if has_pyridine_reactant and has_piperidine_reactant and has_both_in_product:
                        connection_reactions.append((rsmi, depth))
                        print(f"Found connection reaction at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Analyze results to determine if strategy is present
    has_pyridine = len(pyridine_nodes) > 0
    has_piperidine = len(piperidine_nodes) > 0
    has_connection = len(connection_nodes) > 0
    has_connection_reaction = len(connection_reactions) > 0

    # Special case: If the target molecule (depth 0) already contains both rings connected,
    # this is a valid example of the strategy
    if has_connection and any(depth == 0 for _, depth in connection_nodes):
        print("Found linear pyridine-piperidine connection strategy (already connected in target)")
        result = True
        # Structural constraint: co-occurrence_in_molecule(pyridine, piperidine) at last_stage
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "co-occurrence_in_molecule(pyridine, piperidine)",
                "position": "last_stage"
            }
        })

    # Check for linear synthesis pattern:
    # 1. Both scaffolds should be present
    # 2. Either a connection reaction or a connected product should be found
    # 3. The connected structure should appear at a lower depth than individual scaffolds

    if has_pyridine and has_piperidine and (has_connection or has_connection_reaction):
        # Check if the connection appears later in synthesis (lower depth)
        if has_connection:
            min_connection_depth = min(depth for _, depth in connection_nodes)

            # Find scaffolds at greater depths (earlier in synthesis) than the connection
            pyridine_earlier = [
                depth for _, depth in pyridine_nodes if depth > min_connection_depth
            ]
            piperidine_earlier = [
                depth for _, depth in piperidine_nodes if depth > min_connection_depth
            ]

            # If we find individual scaffolds at earlier synthesis stages, this is our pattern
            if pyridine_earlier and piperidine_earlier:
                print("Found linear pyridine-piperidine connection strategy")
                result = True
                # Structural constraint: sequence (pyridine, piperidine) before co-occurrence
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": [
                            "pyridine",
                            "piperidine"
                        ],
                        "after": "co-occurrence_in_molecule(pyridine, piperidine)"
                    }
                })

        # If we found a connection reaction, that's sufficient evidence
        if has_connection_reaction:
            print("Found linear pyridine-piperidine connection strategy via reaction")
            result = True
            # Structural constraint: co-occurrence (pyridine, piperidine, N-arylation)
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "pyridine",
                        "piperidine",
                        "N-arylation"
                    ]
                }
            })

    # Check if we have evidence of separate synthesis of both scaffolds
    # followed by their connection (even if we don't have the exact connection reaction)
    if has_pyridine and has_piperidine and has_connection:
        # Get the minimum depth where we see a connection
        min_connection_depth = min(depth for _, depth in connection_nodes)

        # Check if we have individual scaffolds at greater depths
        isolated_pyridine = any(depth > min_connection_depth for _, depth in pyridine_nodes)
        isolated_piperidine = any(depth > min_connection_depth for _, depth in piperidine_nodes)

        if isolated_pyridine and isolated_piperidine:
            print(
                "Found linear pyridine-piperidine connection strategy (inferred from synthesis path)"
            )
            result = True
            # Structural constraint: sequence (pyridine, piperidine) before co-occurrence
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": [
                        "pyridine",
                        "piperidine"
                    ],
                    "after": "co-occurrence_in_molecule(pyridine, piperidine)"
                }
            })

    return result, findings_json