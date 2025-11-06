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


# Refactoring for Enumeration: Define lists of reactions at the module level
AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Schotten-Baumann_amide",
]

TERMINAL_ALKENE_REACTIONS = [
    "Wittig",
    "Julia Olefination",
    "Heck_terminal_vinyl",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step strategy where a late-stage amide coupling follows an earlier terminal alkene installation. Amide couplings are identified by checking for reactions in the AMIDE_COUPLING_REACTIONS list (e.g., 'Carboxylic acid with primary amine to amide'). Terminal alkene installations are identified by checking for reactions in the TERMINAL_ALKENE_REACTIONS list (e.g., 'Wittig', 'Julia Olefination').
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

    # Initialize tracking variables
    has_amide_coupling = False
    has_terminal_alkene_installation = False

    # Track reaction depths for determining stage of synthesis
    amide_coupling_depth = None
    terminal_alkene_installation_depth = None

    print(f"Starting traversal of route with root SMILES: {route['smiles']}")

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_coupling, has_terminal_alkene_installation
        nonlocal amide_coupling_depth, terminal_alkene_installation_depth, findings_json

        if node["type"] == "mol":
            print(f"Processing molecule node at depth {depth}: {node['smiles'][:30]}...")

        elif node["type"] == "reaction":
            print(f"Processing reaction node at depth {depth}")
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"  Reaction SMILES: {rsmi[:50]}...")

                # Check for amide coupling reaction by reaction type
                for r in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        has_amide_coupling = True
                        amide_coupling_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        print(f"Detected amide coupling at depth {depth} by reaction type")
                        break # Found one, no need to check others

                # Check for terminal alkene installation
                for r in TERMINAL_ALKENE_REACTIONS:
                    if checker.check_reaction(r, rsmi):
                        has_terminal_alkene_installation = True
                        terminal_alkene_installation_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append(r)
                        print(
                            f"Detected terminal alkene installation at depth {depth} by reaction type"
                        )
                        break # Found one, no need to check others

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children recursively
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # If current node is 'mol' (chemical)
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Print summary of findings
    print(
        f"Summary - Amide coupling: {has_amide_coupling}, Terminal alkene: {has_terminal_alkene_installation}"
    )

    result = False
    # Check for the full strategy
    if has_amide_coupling and has_terminal_alkene_installation:
        # Add the co-occurrence constraint if both are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "amide_coupling",
                    "terminal_alkene_installation"
                ],
                "description": "The route must contain at least one reaction from the amide coupling list and one reaction from the terminal alkene installation list."
            }
        })

        # Check relative depths to ensure correct ordering
        if (
            amide_coupling_depth is not None
            and terminal_alkene_installation_depth is not None
        ):

            # Lower depth means later in synthesis (closer to final product)
            if amide_coupling_depth < terminal_alkene_installation_depth:
                print("Detected complete late-stage amide coupling strategy")
                result = True
                # Add the sequence constraint if the order is correct
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "terminal_alkene_installation",
                        "after": "amide_coupling",
                        "description": "A terminal alkene installation must occur in an earlier synthetic step (higher depth) than the amide coupling."
                    }
                })

    return result, findings_json
