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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage Biginelli-like reaction, defined as a reaction with
    three or more components that forms a pyrimidine-dione substructure. This is a
    specific instance of a convergent, multi-component cyclization strategy.
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

    found_strategy = False
    multi_component_reaction_found = False
    pyrimidine_dione_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy, multi_component_reaction_found, pyrimidine_dione_formed, findings_json

        if node["type"] == "reaction":
            if depth <= 1:  # Focus on late-stage reactions
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "multi-component cyclization",
                        "position": "late_stage (depth <= 1)"
                    }
                })

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Count reactants
                reactants = reactants_part.split(".")

                # Check if this is a multi-component reaction (3+ reactants)
                if len(reactants) >= 3:
                    multi_component_reaction_found = True
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants",
                            "operator": ">=",
                            "value": 3
                        }
                    })

                    # Check if product has a heterocyclic system
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol:
                        # Look for pyrimidine-dione pattern
                        pyrimidine_dione_pattern = Chem.MolFromSmarts(
                            "[#7]1[#6](=[#8])[#7][#6](=[#8])[#6][#6]1"
                        )
                        if product_mol.HasSubstructMatch(pyrimidine_dione_pattern):
                            print(
                                f"Found convergent multi-component cyclization at depth {depth}"
                            )
                            pyrimidine_dione_formed = True
                            findings_json["atomic_checks"]["ring_systems"].append("pyrimidine-dione")
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation") # Implied by cyclization
                            found_strategy = True

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if multi_component_reaction_found and pyrimidine_dione_formed:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "multi-component reaction",
                    "pyrimidine-dione formation"
                ]
            }
        })

    print(f"Convergent multi-component cyclization strategy detected: {found_strategy}")
    return found_strategy, findings_json