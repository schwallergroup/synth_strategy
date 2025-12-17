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
    This function detects a convergent synthesis with late-stage amide formation.
    It looks for amide bond formation in the second half of the synthesis
    and checks if multiple fragments are combined.
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

    # Track if we found amide formation and at what depth
    amide_formation_found = False
    amide_formation_depth = -1
    fragment_count = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found, amide_formation_depth, fragment_count, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Check if this is an amide formation reaction
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                acid_chloride_present_flag = False
                amide_in_product_flag = False

                # Check for acid chloride in reactants
                for r in reactants:
                    if r:
                        mol_r = Chem.MolFromSmiles(r)
                        if mol_r is not None and mol_r.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[Cl]")):
                            acid_chloride_present_flag = True
                            if "acid chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("acid chloride")
                            break

                # Check for amide in product
                product_mol = Chem.MolFromSmiles(product) if product else None
                if product_mol is not None and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[N]")
                ):
                    amide_in_product_flag = True
                    if "amide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("amide")

                # If we have acid chloride reactant → amide product, it's an amide formation
                if acid_chloride_present_flag and amide_in_product_flag:
                    amide_formation_found = True
                    amide_formation_depth = depth
                    if "amide_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")
                    print(f"Found amide formation at depth {depth}")

                # Count fragments being combined
                if len(reactants) > 1:
                    fragment_count += len(reactants) - 1
                    print(f"Found {len(reactants)} fragments being combined at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this is a late-stage amide formation in a convergent synthesis
    is_late_stage = amide_formation_depth >= 0 and amide_formation_depth <= max_depth / 2
    is_convergent = fragment_count >= 2

    print(
        f"Amide formation: {amide_formation_found}, Depth: {amide_formation_depth}, Max depth: {max_depth}"
    )
    print(f"Fragment count: {fragment_count}")
    print(f"Is late stage: {is_late_stage}, Is convergent: {is_convergent}")

    result = amide_formation_found and is_late_stage and is_convergent

    # Record structural constraints if conditions are met
    if amide_formation_found and is_late_stage:
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "target": "amide_formation",
                    "position": "late_stage"
                }
            }
        )
    if is_convergent:
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": "fragment_combination_step",
                    "operator": ">=",
                    "value": 2
                }
            }
        )
    if result:
        findings_json["structural_constraints"].append(
            {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "amide_formation",
                        "convergent_synthesis_route"
                    ],
                    "description": "The strategy requires both a late-stage amide formation and a convergent route structure (>= 2 fragment combinations)."
                }
            }
        )

    return result, findings_json
