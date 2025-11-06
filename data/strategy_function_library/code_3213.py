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


HETEROCYCLE_SMARTS = [
    "[n]1[c][n][c][c]1",  # imidazole
    "[s]1[c][n][c][c]1",  # thiazole
    "[o]1[c][c][c][c]1",  # furan
    "[nH]1[c][c][c][c]1", # pyrrole
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage (final step) convergent coupling of two complex fragments (>=10 atoms each),
    where at least one fragment contains a specific heterocycle (imidazole, thiazole, furan, or pyrrole).
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

    found_late_convergent = False
    heterocycle_patterns = [Chem.MolFromSmarts(p) for p in HETEROCYCLE_SMARTS]
    heterocycle_names = ["imidazole", "thiazole", "furan", "pyrrole"]

    def dfs_traverse(node, depth=0):
        nonlocal found_late_convergent, findings_json

        if (
            depth <= 1
            and node["type"] == "reaction"
            and "metadata" in node
            and "mapped_reaction_smiles" in node["metadata"]
        ):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Split reactants
            reactant_smiles = reactants_part.split(".")

            # Only consider reactions with at least 2 reactants
            if len(reactant_smiles) >= 2:
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "reactants_in_coupling_step",
                        "operator": ">=",
                        "value": 2
                    }
                })

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactant_smiles]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and all(r for r in reactant_mols):
                    # Check if at least one reactant has a heterocycle
                    has_heterocycle = False
                    for r_idx, r in enumerate(reactant_mols):
                        if r:
                            for p_idx, p in enumerate(heterocycle_patterns):
                                if r.HasSubstructMatch(p):
                                    has_heterocycle = True
                                    if heterocycle_names[p_idx] not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle_names[p_idx])
                                    break
                            if has_heterocycle: # Found heterocycle in this reactant, no need to check other patterns for this reactant
                                break

                    # Check if reactants are complex (at least 10 atoms)
                    complex_reactants_count = sum(1 for r in reactant_mols if r and r.GetNumAtoms() >= 10)

                    if has_heterocycle and complex_reactants_count >= 2:
                        found_late_convergent = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "convergent_coupling_check",
                                "position": "late_stage"
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "heterocycle_in_reactant",
                                    "at_least_two_complex_reactants"
                                ]
                            }
                        })
                        print(f"Late-stage convergent coupling detected at depth {depth}")
                        print(f"Reactants: {reactant_smiles}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical')
                new_depth = depth + 1
            # If current node is 'reaction', new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    print(f"Late-stage convergent coupling strategy detected: {found_late_convergent}")

    return found_late_convergent, findings_json
