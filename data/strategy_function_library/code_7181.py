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
    Detects if a trifluoromethoxy-containing aryl group is introduced in the final step of the synthesis (depth=1).
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

    trifluoromethoxy_aryl_pattern = Chem.MolFromSmarts("[cX3]:[cX3][OX2][CX4]([FX1])([FX1])[FX1]")

    late_stage_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_introduction, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (final step)
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(trifluoromethoxy_aryl_pattern):
                findings_json["atomic_checks"]["functional_groups"].append("trifluoromethoxy-aryl")
                # Check if trifluoromethoxy aryl was introduced in this reaction
                trifluoromethoxy_in_reactants = False

                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        trifluoromethoxy_aryl_pattern
                    ):
                        trifluoromethoxy_in_reactants = True
                        break

                if not trifluoromethoxy_in_reactants:
                    late_stage_introduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("functional_group_introduction")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "functional_group_introduction",
                            "context": "trifluoromethoxy-aryl",
                            "position": "last_stage"
                        }
                    })
                    print(
                        f"Late-stage introduction of trifluoromethoxy aryl group at depth {depth}"
                    )

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_introduction, findings_json
