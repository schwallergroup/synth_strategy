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
    This function detects SNAr reactions forming diaryl ethers from aryl chlorides.
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

    aryl_chloride_pattern = Chem.MolFromSmarts("[c]-[Cl]")
    diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
    snar_found = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_found, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check if reactant has aryl chloride and product has diaryl ether,
                    # and that the diaryl ether was not present in the reactants.
                    has_aryl_chloride = reactants_mol.HasSubstructMatch(aryl_chloride_pattern)
                    has_product_diaryl_ether = product_mol.HasSubstructMatch(diaryl_ether_pattern)
                    has_reactant_diaryl_ether = reactants_mol.HasSubstructMatch(diaryl_ether_pattern)

                    if has_aryl_chloride:
                        if "aryl chloride" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("aryl chloride")
                    if has_product_diaryl_ether:
                        if "diaryl ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("diaryl ether")

                    if (
                        has_aryl_chloride
                        and has_product_diaryl_ether
                        and not has_reactant_diaryl_ether
                    ):
                        snar_found = True
                        if "SNAr" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("SNAr")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return snar_found, findings_json