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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects reactions that result in a net increase in the number of defined stereocenters.
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

    # Track if we found stereocenter formation
    stereocenter_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal stereocenter_formation_found, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for stereocenter formation
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product_mol = Chem.MolFromSmiles(product_part) if product_part else None

                    if all(reactant_mols) and product_mol:
                        # Count stereocenters in reactants and product
                        reactant_stereocenters = sum(
                            len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
                            for mol in reactant_mols
                        )
                        product_stereocenters = len(
                            Chem.FindMolChiralCenters(product_mol, includeUnassigned=True)
                        )

                        # Check for a net increase in stereocenters in the forward direction.
                        if product_stereocenters > reactant_stereocenters:
                            stereocenter_formation_found = True
                            findings_json["atomic_checks"]["named_reactions"].append("stereocenter_formation")
                except Exception as e:
                    pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Add structural constraint if stereocenter formation was found
    if stereocenter_formation_found:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "stereocenter_formation",
                "operator": ">=",
                "value": 1
            }
        })

    return stereocenter_formation_found, findings_json
