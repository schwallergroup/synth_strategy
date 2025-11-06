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
    Detects if a specific tetrahydroisoquinoline-like scaffold is preserved
    in the product of every reaction step throughout a synthetic route.
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

    # Track scaffold preservation
    scaffold_preserved = True

    # SMARTS pattern for tetrahydroisoquinoline scaffold
    scaffold_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6]2[#6]1[#6][#6][#6][#6]2")

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_preserved, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check if scaffold is preserved
                    if not product_mol.HasSubstructMatch(scaffold_pattern):
                        scaffold_preserved = False
                    else:
                        # If scaffold is preserved in this step, record it as a detected ring system
                        # This is a bit of a conceptual leap, as the strategy is about *non-loss*, 
                        # but for atomic checks, we record what *is* present.
                        # We only add it once to avoid duplicates if it's found in multiple steps.
                        if "tetrahydroisoquinoline" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("tetrahydroisoquinoline")

                except:
                    print("Error processing SMILES in reaction")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if scaffold is preserved throughout the route
    strategy_present = scaffold_preserved

    # If the scaffold was preserved throughout, then the structural constraint is met
    if strategy_present:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_product_without_tetrahydroisoquinoline",
                "operator": "==",
                "value": 0
            }
        })

    print(f"Scaffold preservation strategy detected: {strategy_present}")
    return strategy_present, findings_json