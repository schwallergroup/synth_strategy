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
    This function detects if the synthetic route involves a late-stage
    conversion of a fluoro group to a methoxy group.

    We're looking for a forward reaction where:
    - A reactant has a fluoro group.
    - The product has a methoxy group replacing the fluoro group.
    - The reaction is a late-stage step (depth <= 2).
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

    late_stage_conversion = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_conversion, findings_json
        if late_stage_conversion:  # Stop traversal if already found
            return

        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product_smiles)
                if not product_mol:
                    return

                # Check if product has a methoxy group
                methoxy_patt = Chem.MolFromSmarts("[#6]-[O]-[CH3]")
                if product_mol.HasSubstructMatch(methoxy_patt):
                    if "methoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("methoxy")
                else:
                    return

                # Find a reactant that contains a fluorine atom
                reactant_with_fluoro = None
                for r in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol:
                        for atom in r_mol.GetAtoms():
                            if atom.GetSymbol() == "F":
                                reactant_with_fluoro = r
                                if "fluoro" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("fluoro")
                                break
                        if reactant_with_fluoro:
                            break

                if reactant_with_fluoro:
                    # Verify methoxy-for-fluoro swap by checking atom mappings
                    reactant_mol = Chem.MolFromSmiles(reactant_with_fluoro)
                    methoxy_matches = product_mol.GetSubstructMatches(methoxy_patt)

                    for match in methoxy_matches:
                        # The C attached to the OMe group is at index 0 of the match
                        attachment_atom_prod = product_mol.GetAtomWithIdx(match[0])

                        if attachment_atom_prod.HasProp("molAtomMapNumber"):
                            map_num = attachment_atom_prod.GetProp("molAtomMapNumber")

                            for r_atom in reactant_mol.GetAtoms():
                                if (
                                    r_atom.HasProp("molAtomMapNumber")
                                    and r_atom.GetProp("molAtomMapNumber") == map_num
                                ):
                                    # Check if this mapped atom has a Fluorine neighbor in the reactant
                                    for neighbor in r_atom.GetNeighbors():
                                        if neighbor.GetSymbol() == "F":
                                            late_stage_conversion = True
                                            if "fluoro_to_methoxy_substitution" not in findings_json["atomic_checks"]["named_reactions"]:
                                                findings_json["atomic_checks"]["named_reactions"].append("fluoro_to_methoxy_substitution")
                                            # Add structural constraint if late-stage conversion is detected
                                            findings_json["structural_constraints"].append({
                                                "type": "positional",
                                                "details": {
                                                    "target": "fluoro_to_methoxy_substitution",
                                                    "position": "max_depth_from_product",
                                                    "value": 2
                                                }
                                            })
                                            return  # Found it, stop this branch
                                    break  # Found mapped atom, no need to check others

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_conversion, findings_json
