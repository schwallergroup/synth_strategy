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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis follows a linear strategy
    (as opposed to convergent) by checking if each reaction has only
    one non-reagent reactant.
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

    INHERENTLY_CONVERGENT_REACTIONS = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Heck terminal vinyl",
        "Sonogashira alkyne_aryl halide",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Ugi reaction",
        "A3 coupling",
    ]

    INHERENTLY_LINEAR_REACTIONS = [
        "Wittig reaction",
        "Wittig",
        "Wittig reaction with triphenylphosphorane",
        "Wittig with Phosphonium",
        "Julia Olefination",
        "Grignard with CO2 to carboxylic acid",
        "Grignard from aldehyde to alcohol",
        "Grignard from ketone to alcohol",
    ]

    def is_common_reagent(smiles):
        """
        Identifies common reagents that shouldn't count toward convergence.
        """
        # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for common reagent functional groups
        common_reagent_groups = [
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Magnesium halide",
            "Zinc halide",
            "Alkyl lithium",
            "Tin",
            "Silyl protective group",
            "TMS ether protective group",
            "Silane",
        ]

        for fg in common_reagent_groups:
            if checker.check_fg(fg, smiles):
                if fg not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg)
                print(f"Reagent identified by functional group {fg}: {smiles}")
                return True

        # Common reagent molecules - expanded list
        common_reagents = [
            "CCO",
            "CO",
            "O",
            "[OH-]",
            "[H+]",
            "Cl",
            "Br",
            "I",
            "F",
            "CC(=O)O",
            "CC(=O)[O-]",
            "CN",
            "CS(=O)(=O)O",
            "CS(=O)(=O)[O-]",
            "C[N+](C)(C)C",
            "B(O)O",
            "B(OH)3",
            "P(OCC)(OCC)OCC",
            "P(=O)(OCC)(OCC)OCC",
            "CC(C)O",
            "CCOC(=O)C",
            "CC(=O)C",
            "O=C=O",
            "N#N",
            "N=[N+]=[N-]",
            "[Na+]",
            "[K+]",
            "[Li+]",
            "c1ccccc1",
            "CC(=O)OC(=O)C",
            "ClC(=O)C",
            "BrC(=O)C",
            "IC(=O)C",
            "FC(=O)C",
            "O=S(=O)(O)O",
            "[O-]S(=O)(=O)[O-]",
            "N",
            "NN",
        ]

        # Normalize SMILES for comparison
        norm_smiles = Chem.MolToSmiles(mol)
        for reagent in common_reagents:
            reagent_mol = Chem.MolFromSmiles(reagent)
            if reagent_mol and Chem.MolToSmiles(reagent_mol) == norm_smiles:
                print(f"Common reagent molecule identified: {smiles}")
                return True

        return False

    def is_inherently_convergent_reaction(rsmi):
        """Checks if the reaction type is one of the specified inherently convergent reactions."""
        for rxn_type in INHERENTLY_CONVERGENT_REACTIONS:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                print(f"Inherently convergent reaction detected: {rxn_type}")
                return True

        return False

    def is_inherently_linear_reaction(rsmi):
        """Checks if the reaction type is one of the specified inherently linear reactions."""
        for rxn_type in INHERENTLY_LINEAR_REACTIONS:
            if checker.check_reaction(rxn_type, rsmi):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                print(f"Inherently linear reaction detected despite multiple reactants: {rxn_type}")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, findings_json

        # This print statement is for debugging/demonstration of depth
        # print(f"Node type: {node['type']}, Depth: {depth}")

        if node["type"] == "reaction" and is_linear:  # Skip if already non-linear
            # Extract reactants
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if reaction type is inherently convergent
                if is_inherently_convergent_reaction(rsmi):
                    print(f"Non-linear step detected: inherently convergent reaction")
                    is_linear = False
                    if {"type": "negation", "details": {"target": "presence of an inherently convergent reaction", "description": "The route must not contain any reaction from a predefined list of reactions considered inherently convergent (e.g., Suzuki, Negishi)."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "presence of an inherently convergent reaction", "description": "The route must not contain any reaction from a predefined list of reactions considered inherently convergent (e.g., Suzuki, Negishi)."}})
                    # No return here, continue to traverse children to find more non-linearities

                # Check if reaction type is inherently linear despite multiple reactants
                if is_inherently_linear_reaction(rsmi):
                    print(f"Linear step maintained: reaction is inherently linear")
                    # Continue with traversal
                else:
                    # Count significant reactants (excluding common reagents)
                    significant_reactants = []
                    for reactant in reactants:
                        if not is_common_reagent(reactant):
                            significant_reactants.append(reactant)

                    if len(significant_reactants) > 1:
                        print(
                            f"Non-linear step detected with {len(significant_reactants)} significant reactants"
                        )
                        for r in significant_reactants:
                            print(f"  - {r}")
                        is_linear = False
                        if {"type": "count", "details": {"target": "significant reactants per reaction", "operator": "<=", "value": 1, "condition": "This constraint is ignored if the reaction is on a predefined list of inherently linear reactions (e.g., Wittig, Grignard).", "description": "A significant reactant is any reactant that is not identified as a common, simple reagent."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "significant reactants per reaction", "operator": "<=", "value": 1, "condition": "This constraint is ignored if the reaction is on a predefined list of inherently linear reactions (e.g., Wittig, Grignard).", "description": "A significant reactant is any reactant that is not identified as a common, simple reagent."}})
                        # No return here, continue to traverse children to find more non-linearities
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root with initial depth 0
    dfs_traverse(route, depth=0)

    return is_linear, findings_json
