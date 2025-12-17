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


SUPPORTED_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for a late-stage fragment coupling strategy using a defined list of cross-coupling reactions.
    A positive hit requires the joining of at least two 'complex' fragments in the final two synthetic steps.
    A fragment is considered 'complex' if it has at least 8 atoms and contains a ring or at least two functional groups.
    The specific reactions checked are enumerated in the SUPPORTED_COUPLING_REACTIONS list, including Suzuki, Buchwald-Hartwig, Negishi, Stille, and others.
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

    found_late_coupling = False

    def has_multiple_fg(mol, smiles):
        """Helper function to check if a molecule has multiple functional groups"""
        fg_count = 0
        detected_fgs = []
        for fg in [
            "Aromatic halide", "Boronic acid", "Boronic ester", "Carboxylic acid",
            "Primary amine", "Secondary amine", "Tertiary amine", "Primary alcohol",
            "Secondary alcohol", "Tertiary alcohol", "Alkyne", "Nitrile",
            "Ester", "Ether", "Ketone", "Aldehyde",
        ]:
            if checker.check_fg(fg, smiles):
                fg_count += 1
                detected_fgs.append(fg)
                if fg_count >= 2:
                    # Record functional groups if multiple are found
                    for d_fg in detected_fgs:
                        if d_fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(d_fg)
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "functional_groups_on_reactant",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    return True
        return False

    def is_complex_fragment(mol, smiles):
        """Helper function to determine if a fragment is complex"""
        if not mol:
            return False
        atom_count = mol.GetNumAtoms()
        has_ring = mol.GetRingInfo().NumRings() > 0
        if has_ring:
            # No specific ring system name is checked, but the presence of a ring is a finding
            # This is a placeholder, as the strategy JSON doesn't specify named ring systems for this check
            if "any_ring_system" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("any_ring_system")

        return atom_count >= 8 and (has_ring or has_multiple_fg(mol, smiles))

    def dfs_traverse(node, depth=0):
        nonlocal found_late_coupling, findings_json

        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                if len(reactants_smiles) >= 2:
                    is_coupling_reaction = False
                    for rxn_type in SUPPORTED_COUPLING_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_coupling_reaction = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "any_supported_cross_coupling_reaction",
                                    "position": "last_three_stages"
                                }
                            })
                            break

                    if is_coupling_reaction:
                        complex_fragments_count = 0
                        for r_smiles in reactants_smiles:
                            try:
                                mol = Chem.MolFromSmiles(r_smiles)
                                if mol and is_complex_fragment(mol, r_smiles):
                                    complex_fragments_count += 1
                            except Exception:
                                continue
                        
                        if complex_fragments_count >= 2:
                            found_late_coupling = True
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "complex_reactants",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
            except Exception:
                pass

        for child in node.get("children", []):
            if found_late_coupling:
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_late_coupling, findings_json
