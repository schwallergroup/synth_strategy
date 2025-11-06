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


AMINE_PROTECTION_REACTIONS = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of primary amine",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of primary amines",
    "Phthalimide protection",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a two-step synthetic sequence where a nitrile is first reduced to a primary amine, which is subsequently protected. The amine protection step is identified from a specific list of named reactions: 'Boc amine protection', 'Boc amine protection explicit', 'Boc amine protection with Boc anhydride', 'Boc amine protection (ethyl Boc)', 'Boc amine protection of primary amine', 'Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N', 'Acylation of primary amines', and 'Phthalimide protection'.
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

    # Track transformations and molecules by path
    transformations = {}
    molecules = {}

    def dfs_traverse(node, path=(), depth=0):
        nonlocal transformations, molecules, findings_json
        if node["type"] == "mol":
            # Store molecule SMILES at this path
            molecules[path] = node["smiles"]
            print(f"Found molecule at path {path}: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at path {path}: {rsmi}")

            # Check for nitrile reduction to amine
            if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                transformations[path] = ("nitrile_to_amine", product, reactants)
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitrile to amine")
                print(f"Found nitrile reduction at path {path}: {rsmi}")

            # Check for amine protection reactions
            elif any(
                checker.check_fg("Primary amine", r) for r in reactants
            ) and not checker.check_fg("Primary amine", product):
                # Check for various amine protection reactions
                if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

                for reaction_type in AMINE_PROTECTION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        transformations[path] = ("amine_to_protected", product, reactants)
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        print(f"Found amine protection ({reaction_type}) at path {path}: {rsmi}")
                        break

        # Continue traversal with updated path
        for i, child in enumerate(node.get("children", [])):
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, path + (i,), new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check for the sequence
    nitrile_paths = []
    protection_paths = []

    for path, transform_info in transformations.items():
        transform_type, product, reactants = transform_info
        if transform_type == "nitrile_to_amine":
            nitrile_paths.append((path, product, reactants))
        elif transform_type == "amine_to_protected":
            protection_paths.append((path, product, reactants))

    # Sort paths by depth (longer path = earlier in synthesis)
    nitrile_paths.sort(key=lambda x: len(x[0]), reverse=True)
    protection_paths.sort(key=lambda x: len(x[0]), reverse=True)

    print(f"Nitrile reduction paths: {nitrile_paths}")
    print(f"Amine protection paths: {protection_paths}")

    # Check if we have both transformations
    if not nitrile_paths or not protection_paths:
        print("Missing one or both transformations")
        return False, findings_json

    result = False
    # Check for correct sequence (nitrile reduction happens before protection)
    for nitrile_path, nitrile_product, nitrile_reactants in nitrile_paths:
        for protection_path, protection_product, protection_reactants in protection_paths:
            # In retrosynthesis, shorter path = later in forward synthesis
            # So protection should have a shorter path than nitrile reduction
            if len(protection_path) < len(nitrile_path):
                print(
                    f"Found potential sequence: nitrile reduction at depth {len(nitrile_path)}, protection at depth {len(protection_path)}"
                )
                result = True
                # Add structural constraint finding
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "before": "Reduction of nitrile to amine",
                        "after_any_of": [
                            "Boc amine protection",
                            "Boc amine protection explicit",
                            "Boc amine protection with Boc anhydride",
                            "Boc amine protection (ethyl Boc)",
                            "Boc amine protection of primary amine",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Acylation of primary amines",
                            "Phthalimide protection"
                        ]
                    }
                })
                return result, findings_json

    print("No valid nitrile\u2192amine\u2192protected amine sequence found")
    return result, findings_json
