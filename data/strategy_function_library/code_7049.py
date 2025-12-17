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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Buchwald-Hartwig",
    "N-arylation",
    "Negishi",
    "Stille reaction",
    "Sonogashira",
    "Heck",
    "Ullmann",
    "Ullmann-Goldberg",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where two complex fragments, each with their own synthetic history, are joined in a late-stage coupling reaction. This check specifically identifies the following coupling reaction types: Suzuki, Buchwald-Hartwig, N-arylation, Negishi, Stille reaction, Sonogashira, Heck, Ullmann, Ullmann-Goldberg, Ullmann condensation, Hiyama-Denmark Coupling, Kumada cross-coupling, Aryllithium cross-coupling, and decarboxylative_coupling.
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

    final_coupling_combines_complex_fragments = False

    def is_complex_fragment(smiles):
        """Helper function to determine if a fragment is complex"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Check for rings
            ring_info = mol.GetRingInfo()
            num_rings = len(ring_info.AtomRings())
            num_atoms = mol.GetNumAtoms(onlyExplicit=True)
            num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])

            # A complex fragment should have rings AND either sufficient size OR multiple heteroatoms
            return num_rings >= 1 and (num_atoms >= 8 or num_heteroatoms >= 2)
        return False

    def is_coupling_reaction(rxn_smiles):
        """Check if the reaction is a coupling reaction"""
        for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"Detected coupling reaction: {rxn_type}")
                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True
        return False

    def has_synthetic_history(node):
        """Check if a node has synthetic history (at least one reaction step)"""
        if "children" in node and len(node["children"]) > 0:
            for child in node["children"]:
                if child["type"] == "reaction":
                    return True
                elif child["type"] == "mol" and "in_stock" in child and not child["in_stock"]:
                    if has_synthetic_history(child):
                        return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal final_coupling_combines_complex_fragments, findings_json

        # Check for convergent synthesis at late stages (final, penultimate, or antepenultimate step)
        if node["type"] == "reaction" and depth <= 2:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is a coupling reaction
                if is_coupling_reaction(rsmi):
                    print(f"Found coupling reaction at depth {depth}: {rsmi}")
                    # Record late_stage positional constraint
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "convergent_coupling_reaction",
                            "position": "late_stage",
                            "max_depth": 2
                        }
                    })

                    reactants = rsmi.split(">")[0].split(".")

                    # Check if we have at least 2 reactants
                    if len(reactants) >= 2:
                        # Get reactant child nodes
                        reactant_children = []
                        for child in node.get("children", []):
                            if child["type"] == "mol" and "smiles" in child:
                                reactant_children.append(child)

                        # Count complex fragments with synthetic history
                        complex_fragments_count = 0

                        # Check each reactant for complexity and synthetic history
                        for child in reactant_children:
                            if is_complex_fragment(child["smiles"]):
                                print(f"Found complex fragment: {child['smiles']}")
                                if has_synthetic_history(child):
                                    complex_fragments_count += 1
                                    print(f"Fragment has synthetic history: {child['smiles']}")

                        # True convergent synthesis requires at least 2 complex fragments
                        if complex_fragments_count >= 2:
                            final_coupling_combines_complex_fragments = True
                            print(
                                f"Convergent synthesis confirmed with {complex_fragments_count} complex fragments at depth {depth}"
                            )
                            # Record count constraint
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "complex_reactants_with_synthetic_history",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Convergent synthesis strategy detected: {final_coupling_combines_complex_fragments}")
    return final_coupling_combines_complex_fragments, findings_json
