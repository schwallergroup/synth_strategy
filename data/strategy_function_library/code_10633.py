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
    Detects if a new chiral center (atomic stereocenter) is introduced in the second half of the synthesis.
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

    total_reactions = 0
    stereo_introductions = []
    result = False

    # First pass: count total reactions
    def count_reactions(node):
        nonlocal total_reactions
        if node["type"] == "reaction":
            total_reactions += 1
        for child in node.get("children", []):
            count_reactions(child)

    count_reactions(route)
    print(f"Total reactions in route: {total_reactions}")

    # Second pass: find where stereocenter is introduced
    def find_stereo_introduction(node, depth=0):
        nonlocal stereo_introductions, findings_json
        if node["type"] == "reaction" and node.get("metadata", {}).get("mapped_reaction_smiles"):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Parse molecules with RDKit
            try:
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if None in reactants_mols or product_mol is None:
                    print(f"Warning: Could not parse molecules at depth {depth}")
                    return

                # Count all stereocenters (assigned and unassigned)
                reactants_all_stereo = sum(
                    len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
                    for mol in reactants_mols
                    if mol is not None
                )
                product_all_stereo = len(
                    Chem.FindMolChiralCenters(product_mol, includeUnassigned=True)
                )

                # Count only assigned stereocenters
                reactants_assigned_stereo = sum(
                    len(Chem.FindMolChiralCenters(mol, includeUnassigned=False))
                    for mol in reactants_mols
                    if mol is not None
                )
                product_assigned_stereo = len(
                    Chem.FindMolChiralCenters(product_mol, includeUnassigned=False)
                )

                print(
                    f"Depth {depth}: Reactants stereocenters: all={reactants_all_stereo}, assigned={reactants_assigned_stereo}"
                )
                print(
                    f"Depth {depth}: Product stereocenters: all={product_all_stereo}, assigned={product_assigned_stereo}"
                )

                # Check for introduction of new stereocenters or assignment of existing ones
                if product_assigned_stereo > reactants_assigned_stereo:
                    print(f"Found stereocenter introduction at depth {depth}")
                    stereo_introductions.append(depth)
                    if "chiral_center_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("chiral_center_formation")

            except Exception as e:
                print(f"Error analyzing stereochemistry at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            find_stereo_introduction(child, new_depth)

    find_stereo_introduction(route)
    print(f"Stereocenter introductions at depths: {stereo_introductions}")

    # Check if any stereocenter was introduced in second half of synthesis
    if stereo_introductions:
        # In retrosynthetic trees, lower depth means later in the synthesis
        # We need at least one stereocenter introduction in the second half
        for depth in stereo_introductions:
            if depth < (total_reactions / 2):
                print(f"Late-stage stereocenter introduction detected at depth {depth}")
                result = True
                # Add the structural constraint if the condition is met
                if {"type": "positional", "details": {"target": "chiral_center_formation", "position": "second_half"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "chiral_center_formation", "position": "second_half"}})
                break # Only need one to satisfy the condition

    if not result:
        print("No late-stage stereocenter introduction detected")
    return result, findings_json