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


CONVERGENT_COUPLING_REACTIONS = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann",
    "Kumada",
    "Hiyama-Denmark Coupling",
    "decarboxylative_coupling"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a convergent synthesis strategy
    where two complex fragments are joined in a late-stage coupling reaction.
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

    convergent_synthesis = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Focus on late-stage reactions (low depth in retrosynthetic tree)
            if depth <= 3:
                # Record late-stage constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "convergent_coupling_reaction",
                        "position": "late_stage",
                        "max_depth": 3
                    }
                })

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Only consider reactions with multiple reactants
                if len(reactants) >= 2:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                        # Check if this is a coupling reaction
                        is_coupling = False
                        for rxn_type in CONVERGENT_COUPLING_REACTIONS:
                            if checker.check_reaction(rxn_type, rsmi):
                                is_coupling = True
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                break

                        # Check complexity of reactants
                        complex_fragments = []
                        for i, r_mol in enumerate(reactant_mols):
                            if r_mol is None:
                                continue

                            # Consider fragments with >8 atoms as potentially complex
                            if r_mol.GetNumAtoms() > 8:
                                # Check for ring structures to confirm complexity
                                has_ring = r_mol.GetRingInfo().NumRings() > 0
                                if has_ring:
                                    # Note: The original strategy JSON doesn't list specific ring systems
                                    # to check, so we just note that a ring was found if it contributes
                                    # to complexity. No specific ring system name to add to findings_json["atomic_checks"]["ring_systems"]
                                    pass

                                # Check for functional groups to confirm complexity
                                has_fg = False
                                for fg in [
                                    "Aromatic halide",
                                    "Boronic acid",
                                    "Boronic ester",
                                    "Carboxylic acid",
                                    "Ester",
                                    "Amide",
                                    "Amine",
                                    "Aromatic alcohol",
                                    "Nitrile",
                                    "Alkyne",
                                    "Alkene",
                                ]:
                                    if checker.check_fg(fg, Chem.MolToSmiles(r_mol)):
                                        has_fg = True
                                        findings_json["atomic_checks"]["functional_groups"].append(fg)

                                if has_ring or has_fg:
                                    complex_fragments.append((i, r_mol))

                        # If we have at least 2 complex fragments being joined in a coupling reaction
                        if len(complex_fragments) >= 2 and is_coupling:
                            convergent_synthesis = True
                            # Record count constraint if met
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "complex_reactants_in_coupling_reaction",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
                    except Exception as e:
                        pass

        # Process children with adjusted depth
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return convergent_synthesis, findings_json