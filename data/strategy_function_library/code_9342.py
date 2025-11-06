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
    Detects if the synthesis employs a late-stage fluorination strategy.
    Specifically, checks if a fluorine atom is introduced in the final reactions (depth ≤ 2).
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

    final_reaction_has_fluorination = False

    def dfs_traverse(node, depth=0):
        nonlocal final_reaction_has_fluorination, findings_json

        if (
            node["type"] == "reaction" and depth <= 2
        ):
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a known fluorination reaction
                if checker.check_reaction("Aromatic fluorination", rsmi):
                    final_reaction_has_fluorination = True
                    findings_json["atomic_checks"]["named_reactions"].append("Aromatic fluorination")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "fluorination", "position": {"depth": {"operator": "<=", "value": 2}}}})
                    return
                elif checker.check_reaction("Fluorination", rsmi):
                    final_reaction_has_fluorination = True
                    findings_json["atomic_checks"]["named_reactions"].append("Fluorination")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "fluorination", "position": {"depth": {"operator": "<=", "value": 2}}}})
                    return

                # Create RDKit molecules
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                ]

                if product_mol:
                    # Count fluorine atoms in product
                    product_f_atoms = len(
                        product_mol.GetSubstructMatches(Chem.MolFromSmarts("[F]"))
                    )

                    # Count fluorine atoms in reactants
                    reactants_f_atoms = sum(
                        len(r.GetSubstructMatches(Chem.MolFromSmarts("[F]")))
                        for r in reactant_mols
                        if r is not None
                    )

                    # Check if fluorine atoms increased
                    if product_f_atoms > reactants_f_atoms:
                        final_reaction_has_fluorination = True
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "fluorination", "position": {"depth": {"operator": "<=", "value": 2}}}})
                        return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return final_reaction_has_fluorination, findings_json
