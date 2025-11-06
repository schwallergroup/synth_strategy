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


from rdkit import Chem

BIARYL_COUPLING_REACTIONS = [
    "Suzuki",
    "Stille",
    "Negishi",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Ullmann condensation",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects a late-stage biaryl coupling, a key indicator of a convergent synthetic strategy. The function first checks if the final synthetic step is a named reaction from the predefined list: Suzuki, Stille, Negishi, Kumada cross-coupling, Hiyama-Denmark Coupling, Ullmann condensation, and Aryllithium cross-coupling. As a fallback, it performs a generic structural analysis to identify the formation of any new C-C single bond between aromatic rings originating from two different reactant molecules."""
    
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    biaryl_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_coupling_found, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late-stage reaction (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str)

                if product and all(r for r in reactants) and len(reactants) >= 2:
                    # First check if this is a known biaryl coupling reaction
                    for name in BIARYL_COUPLING_REACTIONS:
                        if checker.check_reaction(name, rsmi):
                            print(f"Found known biaryl coupling reaction at depth {depth}")
                            biaryl_coupling_found = True
                            if name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(name)
                            break
                    
                    if biaryl_coupling_found and depth <= 1: # If a named reaction was found and it's late-stage
                        if {"type": "positional", "details": {"target": "biaryl_coupling", "position": "depth <= 1"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "biaryl_coupling", "position": "depth <= 1"}})

                    if not biaryl_coupling_found: # Only proceed to structural check if not already found by named reaction
                        # If not a known reaction, check for formation of C-C bond between aromatic rings
                        for bond in product.GetBonds():
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                begin_atom = product.GetAtomWithIdx(bond.GetBeginAtomIdx())
                                end_atom = product.GetAtomWithIdx(bond.GetEndAtomIdx())

                                if (
                                    begin_atom.GetSymbol() == "C"
                                    and begin_atom.GetIsAromatic()
                                    and end_atom.GetSymbol() == "C"
                                    and end_atom.GetIsAromatic()
                                    and begin_atom.HasProp("molAtomMapNumber")
                                    and end_atom.HasProp("molAtomMapNumber")
                                ):

                                    begin_map = begin_atom.GetProp("molAtomMapNumber")
                                    end_map = end_atom.GetProp("molAtomMapNumber")

                                    # Find which reactant each mapped atom belongs to
                                    begin_reactant_idx = -1
                                    end_reactant_idx = -1

                                    for r_idx, r in enumerate(reactants):
                                        for a in r.GetAtoms():
                                            if a.HasProp("molAtomMapNumber"):
                                                if a.GetProp("molAtomMapNumber") == begin_map:
                                                    begin_reactant_idx = r_idx
                                                if a.GetProp("molAtomMapNumber") == end_map:
                                                    end_reactant_idx = r_idx

                                    # Check if atoms come from different reactants
                                    if (
                                        begin_reactant_idx != -1
                                        and end_reactant_idx != -1
                                        and begin_reactant_idx != end_reactant_idx
                                    ):
                                        # Verify this bond doesn't exist in reactants (it's a new bond)
                                        bond_in_reactants = False
                                        for r in reactants:
                                            for r_bond in r.GetBonds():
                                                r_begin = r.GetAtomWithIdx(r_bond.GetBeginAtomIdx())
                                                r_end = r.GetAtomWithIdx(r_bond.GetEndAtomIdx())

                                                if (
                                                    r_begin.HasProp("molAtomMapNumber")
                                                    and r_end.HasProp("molAtomMapNumber")
                                                    and r_begin.GetProp("molAtomMapNumber") == begin_map
                                                    and r_end.GetProp("molAtomMapNumber") == end_map
                                                ):
                                                    bond_in_reactants = True
                                                    break
                                            if bond_in_reactants:
                                                break

                                        if not bond_in_reactants:
                                            print(
                                                f"Found biaryl coupling at depth {depth} - new C-C bond between aromatic rings"
                                            )
                                            biaryl_coupling_found = True
                                            if "biaryl_coupling_structural" not in findings_json["atomic_checks"]["named_reactions"]:
                                                findings_json["atomic_checks"]["named_reactions"].append("biaryl_coupling_structural")
                                            
                                            if depth <= 1: # If structural biaryl coupling was found and it's late-stage
                                                if {"type": "positional", "details": {"target": "biaryl_coupling", "position": "depth <= 1"}} not in findings_json["structural_constraints"]:
                                                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "biaryl_coupling", "position": "depth <= 1"}})
                                            break

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return biaryl_coupling_found, findings_json
