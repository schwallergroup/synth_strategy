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


PROTECTION_REACTION_TYPES = [
    "Boc amine protection",
    "Boc amine protection explicit",
    "Boc amine protection with Boc anhydride",
    "Boc amine protection (ethyl Boc)",
    "Boc amine protection of secondary amine",
    "Boc amine protection of primary amine",
    "Esterification of Carboxylic Acids",
]

DEPROTECTION_REACTION_TYPES = [
    "Boc amine deprotection",
    "Boc amine deprotection of guanidine",
    "Boc amine deprotection to NH-NH2",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the preservation of an alpha-amino acid scaffold in all synthesized intermediates of a route. The scaffold is defined by a SMARTS pattern matching an amine (primary, secondary, or tertiary) and a carboxyl derivative (acid, ester, or amide) on the same alpha-carbon. The analysis ignores common protection/deprotection reagents and commercially available starting materials. This strategy uses the reaction lists `PROTECTION_REACTION_TYPES` and `DEPROTECTION_REACTION_TYPES` to help identify relevant reagents.
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

    all_nodes_have_amino_acid = True

    def has_amino_acid_scaffold(smiles):
        """Check if a molecule has an alpha-amino acid scaffold or derivative."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # SMARTS for an alpha-amino acid derivative:
        # [NX3,NX4+;!$(NC=O)]: any amine (primary, sec, tert, not amide)
        # [CX4;!$(C=O)]: attached to a saturated alpha-carbon
        # [CX3](=[OX1])[OX2,NX3]: which is attached to a carboxyl group (acid, ester, or amide)
        alpha_amino_pattern = "[NX3,NX4+;!$(NC=O)][CX4;!$(C=O)][CX3](=[OX1])[OX2,NX3]"
        if mol.HasSubstructMatch(Chem.MolFromSmarts(alpha_amino_pattern)):
            if "alpha-amino acid scaffold" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("alpha-amino acid scaffold")
            return True
        return False

    def is_protection_deprotection_reagent(smiles):
        """Check if a molecule is a common protection/deprotection reagent."""
        if smiles == "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C":  # Boc anhydride
            return True
        if smiles == "CC(C)(C)OC(=O)Cl":  # Boc-Cl
            return True
        if checker.check_fg("Boc", smiles) and not has_amino_acid_scaffold(smiles):
            if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Boc")
            return True
        return False

    def is_protection_deprotection_reaction(rxn_smiles):
        """Check if a reaction is a protection or deprotection reaction."""
        for rxn_type in PROTECTION_REACTION_TYPES + DEPROTECTION_REACTION_TYPES:
            if checker.check_reaction(rxn_type, rxn_smiles):
                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal all_nodes_have_amino_acid, findings_json

        if node["type"] == "reaction":
            if is_protection_deprotection_reaction(node["smiles"]):
                # The check_reaction call inside is_protection_deprotection_reaction already adds to findings_json
                pass
            for child in node.get("children", []):
                # Depth remains the same when traversing from a reaction node
                dfs_traverse(child, depth)
            return

        if node["type"] == "mol" and "smiles" in node:
            if node.get("in_stock", False):
                return

            if is_protection_deprotection_reagent(node["smiles"]):
                # The check_fg call inside is_protection_deprotection_reagent already adds to findings_json
                return

            if not has_amino_acid_scaffold(node["smiles"]):
                all_nodes_have_amino_acid = False
                # Record the structural constraint violation
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "event": "A synthesized intermediate (excluding stock chemicals and specific reagents) lacks an alpha-amino acid scaffold."
                    }
                })

        for child in node.get("children", []):
            # Depth increases when traversing from a chemical node (or other non-reaction node)
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return all_nodes_have_amino_acid, findings_json
