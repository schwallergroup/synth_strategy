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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects a single reaction step where a nitrile-containing reactant is converted into a lactone-containing product.
    The transformation is confirmed by verifying that the original nitrile carbon becomes part of a ring structure in the product.
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

    nitrile_to_lactone = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_lactone, findings_json
        if nitrile_to_lactone:
            return

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smi = rsmi.split(">")[0].split(".")
                product_smi = rsmi.split(">")[-1]

                try:
                    nitrile_found_in_reactants = False
                    for r in reactants_smi:
                        if checker.check_fg("Nitrile", r):
                            nitrile_found_in_reactants = True
                            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                            break

                    lactone_found_in_product = checker.check_fg("Lactone", product_smi)
                    if lactone_found_in_product:
                        if "Lactone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Lactone")

                    if nitrile_found_in_reactants and lactone_found_in_product:
                        nitrile_atom_maps = set()
                        for reactant_smi in reactants_smi:
                            if checker.check_fg("Nitrile", reactant_smi):
                                reactant_mol = Chem.MolFromSmiles(reactant_smi)
                                if not reactant_mol: continue
                                for atom in reactant_mol.GetAtoms():
                                    if atom.GetSymbol() == 'C' and not atom.GetIsAromatic():
                                        for bond in atom.GetBonds():
                                            other_atom = bond.GetOtherAtom(atom)
                                            if other_atom.GetSymbol() == 'N' and bond.GetBondType() == Chem.BondType.TRIPLE:
                                                if atom.HasProp('molAtomMapNumber'):
                                                    nitrile_atom_maps.add(atom.GetProp('molAtomMapNumber'))
                                                break

                        if not nitrile_atom_maps:
                            return

                        product_mol = Chem.MolFromSmiles(product_smi)
                        if not product_mol: return
                        
                        ring_carbon_maps = set()
                        for atom in product_mol.GetAtoms():
                            if atom.IsInRing() and atom.GetSymbol() == 'C' and atom.HasProp('molAtomMapNumber'):
                                ring_carbon_maps.add(atom.GetProp('molAtomMapNumber'))

                        if nitrile_atom_maps.intersection(ring_carbon_maps):
                            nitrile_to_lactone = True
                            # Add structural constraint if the condition is met
                            structural_constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Nitrile",
                                        "Lactone"
                                    ],
                                    "scope": "reaction_step",
                                    "relationship": "reactant_to_product"
                                }
                            }
                            if structural_constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_obj)

                except Exception:
                    pass

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return nitrile_to_lactone, findings_json