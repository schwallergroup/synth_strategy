from typing import Tuple, Dict, List
import copy
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


from rdkit import Chem

# A list of common coupling reaction names for identification.
COUPLING_REACTION_NAMES = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann",
    "Kumada",
    "Hiyama-Denmark",
]

# A list of SMILES for common reagents that should not be counted as significant reactants.
COMMON_REAGENT_SMILES = [
    # Solvents, Acids, Bases, and simple reagents
    "O=C(O)C", "CC(=O)O", "O", "CO", "CCO", "[H][H]", "C(=O)=O", "N#N",
    "ClC(=O)C", "CC(=O)Cl", "CC(C)(C)OC(=O)O", "CN(C)C=O", "CS(=O)(=O)O",
    "CC(=O)OC(=O)C", "O=S(=O)(O)O", "Cl", "Br", "I", "F", "N",
    "CC(=O)OC(C)(C)C", "CC(=O)O[C]([CH3])=[O]", "CN1CCCC1=O", "CS(=O)(=O)Cl",
    "CC(C)(C)OC(=O)Cl", "CC(C)C[N+](C)(C)C", "C[N+](C)(C)C", "O=C1CCC(=O)N1",
    "O=C(Cl)Cl", "O=C(O)C(F)(F)F", "O=C(O)C=O", "O=S(=O)(Cl)Cl", "O=P(Cl)(Cl)Cl",
'N#CC(=O)OCC', "O=C(O)CO", "O=C(O)CCl", "O=C(O)CBr", "O=C(O)CF", "O=C(O)CI",
    "CC#N", "C1CCOC1", "C1CCOCC1", "CC(C)=O", "COCCOC", "CC(C)O", "[OH-]",
    "[H+]", "[Na+]", "[K+]", "C[Si](C)(C)Cl", "C[Si](C)(C)OC(C)(C)C",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis is non-linear (convergent) by identifying fragment coupling steps.
    A step is considered a fragment coupling if it is a known coupling reaction (from the COUPLING_REACTION_NAMES list)
    or if it involves more than one significant reactant. Reactants are considered significant if they are not
    found in the COMMON_REAGENT_SMILES list or classified as simple acylating agents.
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

    is_linear = True

    def is_reagent(smiles):
        """Helper function to determine if a molecule is likely a reagent"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return True

        # Check if it's a common reagent by SMILES
        if smiles in COMMON_REAGENT_SMILES:
            return True

        # Check for small acylating agents by pattern
        if checker.check_fg("Acyl halide", smiles) and mol.GetNumHeavyAtoms() <= 8:
            if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
            return True

        if checker.check_fg("Anhydride", smiles) and mol.GetNumHeavyAtoms() <= 10:
            if "Anhydride" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Anhydride")
            return True

        return False

    def is_coupling_reaction(rsmi):
        """Helper function to identify if a reaction is a coupling reaction"""
        for rxn_type in COUPLING_REACTION_NAMES:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, findings_json

        if not is_linear:  # Early return if we already know it's not linear
            return

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a known coupling reaction
            if is_coupling_reaction(rsmi):
                is_linear = False
                # Record the structural constraint for negation
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "targets": [
                            "Suzuki", "Negishi", "Stille", "Heck", "Sonogashira",
                            "Buchwald-Hartwig", "Ullmann", "Kumada", "Hiyama-Denmark"
                        ],
                        "description": "The route must not contain any of the specified named coupling reactions to be considered linear."
                    }
                })
                return

            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Filter out molecules that are likely reagents
            significant_reactants = [r for r in reactants if not is_reagent(r)]

            if len(significant_reactants) > 1:
                is_linear = False
                # Record the structural constraint for count
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "significant_reactants_per_step",
                        "operator": "<=",
                        "value": 1,
                        "description": "Each reaction step must have at most one significant reactant. Reactants are considered insignificant if they are common reagents or simple acylating agents (Acyl halide, Anhydride with size limits)."
                    }
                })

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., chemical)
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    return is_linear, findings_json