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


NITROGEN_CONTAINING_HETEROCYCLES = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "imidazole", "oxazole", "thiazole", "pyrrole", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "benzimidazole",
    "benzoxazole", "benzothiazole", "piperidine", "piperazine", "morpholine",
    "thiomorpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage convergent synthesis (final two steps) via N-alkylation.
    The strategy requires the coupling of two complex fragments (defined as having >15 atoms or >1 ring).
    One fragment must be a secondary amine, and the other must be a halo-substituted nitrogen-containing heterocycle.
    The list of recognized heterocycles is defined in `NITROGEN_CONTAINING_HETEROCYCLES`.
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

    found_pattern = False

    def is_complex_fragment(smiles):
        """Check if a molecule is a complex fragment (has multiple rings or >15 atoms)"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            ring_info = mol.GetRingInfo()
            ring_count = ring_info.NumRings()
            atom_count = mol.GetNumAtoms()
            return ring_count > 1 or atom_count > 15
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, findings_json

        if node["type"] == "reaction" and depth <= 2:
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) >= 2:
                    is_n_alkylation = checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )

                    if is_n_alkylation:
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of secondary amines with alkyl halides")
                        # Positional constraint: late_stage
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "N-alkylation of secondary amines with alkyl halides",
                                "position": "late_stage"
                            }
                        })

                    if not is_n_alkylation:
                        return

                    amine_reactant = None
                    halide_heterocycle_reactant = None

                    for reactant in reactants:
                        if checker.check_fg("Secondary amine", reactant):
                            amine_reactant = reactant
                            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                        has_halide = False
                        for halide_type in ["Primary halide", "Secondary halide", "Tertiary halide", "Aromatic halide"]:
                            if checker.check_fg(halide_type, reactant):
                                has_halide = True
                                if halide_type not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(halide_type)

                        has_heterocycle = False
                        for ring in NITROGEN_CONTAINING_HETEROCYCLES:
                            if checker.check_ring(ring, reactant):
                                has_heterocycle = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)

                        if has_halide and has_heterocycle:
                            halide_heterocycle_reactant = reactant

                    if amine_reactant and halide_heterocycle_reactant:
                        amine_complex = is_complex_fragment(amine_reactant)
                        heterocycle_complex = is_complex_fragment(halide_heterocycle_reactant)

                        if amine_complex and heterocycle_complex:
                            found_pattern = True
                            # Co-occurrence constraint
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "reactant_with_secondary_amine",
                                        "reactant_with_halide_and_nitrogen_heterocycle"
                                    ]
                                }
                            })

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_pattern, findings_json
