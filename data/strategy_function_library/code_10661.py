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


CN_COUPLING_REACTIONS = [
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Goldberg coupling",
    "Ullmann-Goldberg Substitution amine",
    "N-arylation_heterocycles",
    "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent C-N coupling step that joins two complex fragments. The reaction must be one of a specific list of named C-N coupling reactions (e.g., Buchwald-Hartwig, Ullmann, Goldberg), and both reacting fragments must be sufficiently complex (containing a ring or >7 heavy atoms).
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

    convergent_cn_coupling = False

    # Define the structural constraint object to be added when met
    structural_constraint_cn_coupling = {
      "type": "co-occurrence",
      "details": {
        "description": "A single reaction step must be a C-N coupling that joins a complex nitrogen-containing reactant with a complex halide-containing reactant.",
        "targets": [
          "C-N coupling reaction",
          "complex_nitrogen_reactant",
          "complex_halide_reactant"
        ]
      }
    }

    def dfs_traverse(node, depth=0):
        nonlocal convergent_cn_coupling, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if len(reactants_smiles) >= 2:
                    is_cn_coupling = False
                    for reaction_type in CN_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_cn_coupling = True
                            if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            break

                    if is_cn_coupling:
                        reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product = Chem.MolFromSmiles(product_smiles)

                        if product is not None:
                            complex_reactants = []
                            for i, r in enumerate(reactants):
                                if r is not None:
                                    has_ring = len(Chem.GetSSSR(r)) > 0
                                    if has_ring:
                                        # Note: RDKit's GetSSSR returns a tuple of rings, not their names.
                                        # We can't easily get specific ring system names here without more complex logic.
                                        # For now, we'll just note the presence of a ring if it's part of the complexity check.
                                        pass # No specific ring system name to add to findings_json["atomic_checks"]["ring_systems"]
                                    if has_ring or r.GetNumAtoms() > 7:
                                        complex_reactants.append(i)

                            if len(complex_reactants) >= 2:
                                n_reactants = []
                                halide_reactants = []

                                for i, r_idx in enumerate(reactants):
                                    if r_idx is not None:
                                        r_smiles = reactants_smiles[i]
                                        fg_names = [
                                            "Primary amine",
                                            "Secondary amine",
                                            "Tertiary amine",
                                            "Aniline",
                                            "Pyrrole",
                                            "Pyridine",
                                            "Imidazole",
                                        ]
                                        for fg_name in fg_names:
                                            if checker.check_fg(fg_name, r_smiles):
                                                n_reactants.append(i)
                                                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                                break # Found an N-containing FG, move to next reactant

                                        if checker.check_fg("Aromatic halide", r_smiles):
                                            halide_reactants.append(i)
                                            if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                                findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                                if n_reactants and halide_reactants:
                                    for n_idx in n_reactants:
                                        for h_idx in halide_reactants:
                                            if (
                                                n_idx != h_idx
                                                and n_idx in complex_reactants
                                                and h_idx in complex_reactants
                                            ):
                                                convergent_cn_coupling = True
                                                # Add the structural constraint if it's not already there
                                                if structural_constraint_cn_coupling not in findings_json["structural_constraints"]:
                                                    findings_json["structural_constraints"].append(structural_constraint_cn_coupling)
                                                break
                                        if convergent_cn_coupling:
                                            break
            except Exception:
                pass # traceback.print_exc() # For debugging

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return convergent_cn_coupling, findings_json
