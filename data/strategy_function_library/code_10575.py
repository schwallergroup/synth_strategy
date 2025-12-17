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


# Refactored module-level constant
CONVERGENT_REACTION_TYPES = [
    # Cross-coupling reactions
    "Suzuki", "Heck", "Sonogashira", "Stille", "Negishi",
    "Buchwald-Hartwig", "N-arylation", "Ullmann", "Ullmann-Goldberg",
    "Goldberg", "Chan-Lam", "Kumada", "Hiyama-Denmark",
    "Aryllithium cross-coupling",
    # Other multi-component/convergent reactions
    "Diels-Alder", "Aldol condensation", "Wittig", "Ugi", "A3 coupling",
    "Petasis", "Michael addition", "Reductive amination",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects if a synthesis is non-linear (i.e., convergent). A route is flagged as non-linear if any reaction step involves more than one 'significant' reactant. A 'significant' reactant is determined by heuristics based on heavy atom count and structure, aiming to exclude common reagents. An exception is made for known convergent reaction types (e.g., Suzuki, Diels-Alder), which are permitted to have up to two significant reactants. The specific reactions checked are defined in the CONVERGENT_REACTION_TYPES list."""
    is_linear = True

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Common reagents and small molecules that shouldn't count toward convergence
    common_reagents = [
        "O", "O=O", "[OH-]", "[H+]", "Cl", "Br", "I", "F", "[Na+]", "[K+]",
        "CC(=O)O", "CCO", "CO", "CC(C)O", "CN", "CS", "C[NH2]", "CC#N",
        "CC(=O)[O-]", "O=C=O", "N#N", "[H][H]", "C1CCOC1", "C1CCCCC1",
        "ClCCl", "BrCBr", "c1ccccc1", "CC(=O)OC(=O)C", "B(O)O", "B(OH)O",
        "B(OH)2", "CC(=O)Cl", "C(=O)Cl", "S(=O)(=O)Cl", "P(Cl)Cl",
        "P(=O)(Cl)Cl", "CN(C)C", "CC=O", "C=O", "C[MgBr]", "C[MgCl]",
        "C[Li]", "C[ZnBr]", "C[ZnCl]", "OB(O)O", "OB(OH)2", "CCN", "CCNCC",
        "CC(C)N", "CC(N)C", "CS(=O)(=O)O", "CS(=O)(=O)[O-]", "CC(O)=O",
        "CC(=O)OC", "CC(=O)OCC",
    ]

    def is_significant_reactant(smiles):
        """Determine if a reactant is significant for convergence analysis"""
        if smiles in common_reagents:
            return False
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False
            heavy_atom_count = mol.GetNumHeavyAtoms()
            if (
                checker.check_fg("Magnesium halide", smiles)
                or checker.check_fg("Alkyl lithium", smiles)
                or "Sn" in smiles
                or "Zn" in smiles
            ):
                if checker.check_fg("Magnesium halide", smiles) and "Magnesium halide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Magnesium halide")
                if checker.check_fg("Alkyl lithium", smiles) and "Alkyl lithium" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Alkyl lithium")
                if heavy_atom_count <= 8:
                    return False
            if checker.check_fg("Boronic acid", smiles) or checker.check_fg(
                "Boronic ester", smiles
            ):
                if checker.check_fg("Boronic acid", smiles) and "Boronic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
                if checker.check_fg("Boronic ester", smiles) and "Boronic ester" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boronic ester")
                if heavy_atom_count <= 12:
                    return False
            return heavy_atom_count > 5
        except Exception:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, findings_json

        if node["type"] == "reaction" and is_linear:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                significant_reactants = [
                    r for r in reactants if is_significant_reactant(r)
                ]

                is_convergent_type = False
                for rxn_type in CONVERGENT_REACTION_TYPES:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_convergent_type = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                max_allowed = 2 if is_convergent_type else 1

                if len(significant_reactants) > max_allowed:
                    is_linear = False

                # Record structural constraints based on the outcome
                if is_convergent_type:
                    if len(significant_reactants) <= 2:
                        # This constraint was met for a convergent reaction
                        constraint_obj = {"type": "count", "details": {"target": "significant_reactants_per_convergent_reaction_step", "operator": "<=", "value": 2}}
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)
                else:
                    if len(significant_reactants) <= 1:
                        # This constraint was met for a non-convergent reaction
                        constraint_obj = {"type": "count", "details": {"target": "significant_reactants_per_non_convergent_reaction_step", "operator": "<=", "value": 1}}
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)

            except Exception:
                # Silently fail on malformed data, assuming linearity as a safe default.
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return is_linear, findings_json