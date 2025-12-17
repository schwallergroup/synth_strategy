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


from rdkit import Chem

RINGS_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "pyrrole", "pyridine",
    "pyrazole", "imidazole", "oxazole", "thiazole", "pyrimidine", "pyrazine",
    "pyridazine", "triazole", "tetrazole", "pyrrolidine", "piperidine",
    "piperazine", "morpholine", "thiomorpholine", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "thiophene",
    "thiopyran", "benzothiophene", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "cyclopropane", "cyclobutane", "cyclopentane",
    "cyclohexane", "cycloheptane", "cyclooctane", "benzene", "naphthalene",
    "anthracene", "benzoxazole", "benzothiazole", "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves introducing unsaturation into a ring system.
    It looks for reactions where a saturated ring is converted to one containing an alkene.

    Since we're traversing the route retrosynthetically, we're actually looking for reactions
    where an unsaturated ring in the product becomes more saturated in the reactants.
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

    found_transformation = False

    def count_ring_double_bonds(mol):
        """Helper function to count double bonds in rings of a molecule"""
        if not mol:
            return 0

        ring_double_bonds = 0
        ring_info = mol.GetRingInfo()
        ring_atoms_sets = [set(ring) for ring in ring_info.AtomRings()]

        if not ring_atoms_sets:
            return 0

        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # Check if this bond is in any ring
                for ring_atoms in ring_atoms_sets:
                    if begin_idx in ring_atoms and end_idx in ring_atoms:
                        ring_double_bonds += 1
                        break

        return ring_double_bonds

    def dfs_traverse(node, depth=0):
        nonlocal found_transformation, findings_json

        if found_transformation:
            return  # Early exit if we already found what we're looking for

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            try:
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if not product_mol or not all(reactants_mols):
                    print("Could not parse all molecules in the reaction")
                    return

                # Check if any of the molecules contain rings from our list
                has_ring_of_interest = False
                rings_found_in_step = set()
                for ring_type in RINGS_OF_INTEREST:
                    if checker.check_ring(ring_type, product_smiles):
                        has_ring_of_interest = True
                        rings_found_in_step.add(ring_type)
                    for r_smi in reactants_smiles.split("."):
                        if checker.check_ring(ring_type, r_smi):
                            has_ring_of_interest = True
                            rings_found_in_step.add(ring_type)

                if not has_ring_of_interest:
                    print("No rings of interest found in this reaction")
                    return

                # Add found rings to findings_json
                for r_type in rings_found_in_step:
                    if r_type not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(r_type)

                # Count double bonds in rings
                product_ring_double_bonds = count_ring_double_bonds(product_mol)

                reactant_ring_double_bonds = sum(
                    count_ring_double_bonds(mol) for mol in reactants_mols if mol
                )

                print(f"Product ring double bonds: {product_ring_double_bonds}")
                print(f"Reactant ring double bonds: {reactant_ring_double_bonds}")

                # A forward reaction introduces unsaturation if the product has more ring double bonds.
                if product_ring_double_bonds > reactant_ring_double_bonds:
                    found_transformation = True
                    print(f"Confirmed ring unsaturation introduction: {rsmi}")
                    if "ring_unsaturation_on_specified_ring" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_unsaturation_on_specified_ring")
                    return

            except Exception as e:
                print(f"Error analyzing ring unsaturation: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    if found_transformation:
        # Add the structural constraint if the transformation was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ring_unsaturation_on_specified_ring",
                "operator": ">=",
                "value": 1
            }
        })

    return found_transformation, findings_json
