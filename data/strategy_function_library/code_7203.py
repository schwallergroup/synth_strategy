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


DEPROTECTION_REACTIONS = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Tert-butyl deprotection of amine",
    "Alcohol deprotection from silyl ethers",
    "Alcohol deprotection from silyl ethers (double)",
    "Alcohol deprotection from silyl ethers (diol)",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydroxyl benzyl deprotection",
    "Carboxyl benzyl deprotection",
    "Cleavage of methoxy ethers to alcohols",
    "Cleavage of alkoxy ethers to alcohols",
    "COOH ethyl deprotection",
    "N-glutarimide deprotection",
    "Phthalimide deprotection",
    "TMS deprotection from alkyne",
    "Ether cleavage to primary alcohol",
]

AROMATIC_FUNCTIONALIZATIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
    "Aromatic nitration with HNO3",
    "Aromatic nitration with NO3 salt",
    "Aromatic nitration with NO2 salt",
    "Aromatic nitration with alkyl NO2",
    "Friedel-Crafts acylation",
    "Friedel-Crafts alkylation",
    "Friedel-Crafts alkylation with halide",
    "Aromatic hydroxylation",
    "Aromatic sulfonyl chlorination",
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Heck terminal vinyl",
    "Oxidative Heck reaction",
    "Oxidative Heck reaction with vinyl ester",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Minisci (para)",
    "Minisci (ortho)",
    "Minisci (para-cyanide)",
    "Minisci (ortho-cyanide)",
    "Minisci-like halide substitution",
    "Catellani reaction ortho",
    "Catellani reaction para",
    "arylhydroxylation",
    "Chan-Lam alcohol",
    "Chan-Lam amine",
    "Chan-Lam etherification",
    "Williamson Ether Synthesis",
    "Ullmann condensation",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy that combines sequential aromatic
    functionalization with a late-stage deprotection. The reaction types for both
    categories are identified from the module-level constants `AROMATIC_FUNCTIONALIZATIONS`
    and `DEPROTECTION_REACTIONS`.
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

    # Track deprotection reactions with their depths
    deprotection_depths = []
    # Track aromatic functionalization reactions with their depths
    aromatic_func_depths = []
    # Track total depth of the route
    max_depth = [0]

    def is_aromatic_mol(mol_smiles):
        """Check if a molecule contains aromatic rings using RDKit"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                for atom in mol.GetAtoms():
                    if atom.GetIsAromatic():
                        return True
            return False
        except:
            return False

    def is_deprotection_reaction(rsmi):
        """Check if a reaction is a deprotection reaction"""
        for rxn_type in DEPROTECTION_REACTIONS:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"DEBUG: Identified deprotection reaction: {rxn_type}")
                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def is_aromatic_functionalization(rsmi):
        """Check if a reaction is an aromatic functionalization"""
        for rxn_type in AROMATIC_FUNCTIONALIZATIONS:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"DEBUG: Identified aromatic functionalization: {rxn_type}")
                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def traverse_route(node, depth=0):
        # Update max depth
        max_depth[0] = max(max_depth[0], depth)

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"DEBUG: Analyzing reaction at depth {depth}: {rsmi}")

                # Check for deprotection reactions
                if is_deprotection_reaction(rsmi):
                    deprotection_depths.append(depth)
                    print(f"DEBUG: Found deprotection reaction at depth {depth}")

                # Check for aromatic functionalization reactions
                if is_aromatic_functionalization(rsmi):
                    aromatic_func_depths.append(depth)
                    print(f"DEBUG: Found aromatic functionalization at depth {depth}")

            except Exception as e:
                print(f"DEBUG: Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                traverse_route(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                traverse_route(child, depth + 1)

    # Start traversal from the root
    traverse_route(route)

    print(f"DEBUG: Deprotection depths: {deprotection_depths}")
    print(f"DEBUG: Aromatic functionalization depths: {aromatic_func_depths}")
    print(f"DEBUG: Max depth of route: {max_depth[0]}")

    # Check if we have both strategies
    has_late_deprotection = False
    has_sequential_aromatic_functionalization = False
    result = False

    # Determine late-stage threshold based on route depth
    late_stage_threshold = min(3, max(1, max_depth[0] // 3))
    print(f"DEBUG: Late-stage threshold: {late_stage_threshold}")

    # Check for late-stage deprotection (lower depth values)
    if deprotection_depths:
        min_deprotection_depth = min(deprotection_depths)
        # Consider it late-stage if it's in the first few steps
        if min_deprotection_depth <= late_stage_threshold:
            has_late_deprotection = True
            print(f"DEBUG: Found late-stage deprotection at depth {min_deprotection_depth}")
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "deprotection",
                    "position": "late_stage"
                }
            })

    # Check for sequential aromatic functionalization
    # Need at least 2 aromatic functionalization reactions
    if len(aromatic_func_depths) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "aromatic_functionalization",
                "operator": ">=",
                "value": 2
            }
        })
        # Sort depths to check if they're sequential
        sorted_depths = sorted(aromatic_func_depths)
        print(f"DEBUG: Sorted aromatic functionalization depths: {sorted_depths}")

        # Check if any two functionalization steps are sequential (within 2 steps of each other)
        for i in range(len(sorted_depths) - 1):
            if sorted_depths[i + 1] - sorted_depths[i] <= 2:  # Allow 1 step in between
                has_sequential_aromatic_functionalization = True
                print(
                    f"DEBUG: Found sequential aromatic functionalizations at depths {sorted_depths[i]} and {sorted_depths[i+1]}"
                )
                break

    # Check if deprotection occurs after (or at same depth as) aromatic functionalizations in retrosynthesis
    # This means deprotection is at a lower or equal depth number
    if has_late_deprotection and has_sequential_aromatic_functionalization:
        min_deprotection = min(deprotection_depths)

        print(
            f"DEBUG: Min deprotection depth: {min_deprotection}, Aromatic depths: {aromatic_func_depths}"
        )

        # In retrosynthesis, we want the deprotection to be at a lower or equal depth
        # compared to at least one aromatic functionalization
        if min_deprotection <= max(aromatic_func_depths):
            print(
                "DEBUG: Detected combined strategy: sequential aromatic functionalization with late-stage deprotection"
            )
            result = True
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "deprotection",
                        "aromatic_functionalization"
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "aromatic_functionalization",
                    "after": "deprotection"
                }
            })
        else:
            print("DEBUG: Deprotection occurs too early in the synthesis")
    elif not has_late_deprotection:
        print("DEBUG: No late-stage deprotection found")
    elif not has_sequential_aromatic_functionalization:
        print("DEBUG: No sequential aromatic functionalization found")

    return result, findings_json
