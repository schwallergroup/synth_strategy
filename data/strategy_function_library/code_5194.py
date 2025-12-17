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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves the introduction of a
    trifluoromethyl-containing aromatic fragment.
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

    trifluoromethyl_pattern = Chem.MolFromSmarts("c-[C]([F])([F])[F]")
    trifluoromethyl_introduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_introduction_detected, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        product_mol
                        and len(product_mol.GetSubstructMatches(trifluoromethyl_pattern)) > 0
                    ):
                        # Check if any reactant contains trifluoromethyl group
                        reactant_has_trifluoromethyl = False

                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if (
                                reactant_mol
                                and len(reactant_mol.GetSubstructMatches(trifluoromethyl_pattern))
                                > 0
                            ):
                                reactant_has_trifluoromethyl = True
                                break

                        # If product has trifluoromethyl but no reactant does, it's introduced
                        if not reactant_has_trifluoromethyl:
                            trifluoromethyl_introduction_detected = True
                            findings_json["atomic_checks"]["functional_groups"].append("aromatic trifluoromethyl")
                            findings_json["atomic_checks"]["named_reactions"].append("functional_group_introduction")
                except:
                    print("Error processing SMILES in trifluoromethyl detection")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return trifluoromethyl_introduction_detected, findings_json