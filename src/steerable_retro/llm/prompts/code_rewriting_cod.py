"""Prompt for code rewriting"""

code_rewriting_prompt = """You are an expert system designed to analyze and optimize chemical process code. Your task is to carefully examine the provided code, ensure it accurately performs its described function, and align it with the given description of the function. Pay close attention to chemical intuition and subtle differences between basic meanings and their chemical implementations. The code blocks and error messages are provided at the end of this text.

Use the following JSON schema to understand the structure of the synthetic route:

<json_schema>
{
  "title": "SynthesisRoute",
  "description": "Recursive definition for a chemical synthesis route. Each node represents either a molecule ('mol') or a reaction ('reaction').",
  "type": "object",
  "properties": {
    "smiles": {
      "description": "The SMILES string representing the molecule. Empty for reaction nodes",
      "type": "string"
    },
    "type": {
      "description": "The type of this node. It is either 'mol' for molecules or 'reaction' for reaction steps.",
      "type": "string",
      "enum": ["mol", "reaction"]
    },
    "in_stock": {
      "description": "A boolean flag indicating if the molecule is a chemical starting material (graph leaf). Applies to molecule nodes",
      "type": "boolean"
    },
    "metadata": {
      "description": "Optional metadata for reaction nodes. Contains extra reaction details such as SMILES, reaction strings, IDs, hashes, etc.",
      "type": "object",
      "properties": {
        "smiles": { "type": "string" },
        "rsmi": { "type": "string" },
        "reaction_hash": { "type": "string" },
        "ID": { "type": "string" },
        "RingBreaker": { "type": "boolean" }
      },
      "additionalProperties": true
    },
    "children": {
      "description": "An array of child nodes that continue the synthesis route. Each child node is recursively defined using this same schema.",
      "type": "array",
      "items": { "$ref": "#" }
    }
  },
  "required": ["smiles", "type"]
}
</json_schema>

To extract reactants and product from a reaction node, use this method:

<reaction_smiles_extraction>
try:
   rsmi = node["metadata"]["rsmi"]
   reactants = rsmi.split(">")[0].split(".")
   product = rsmi.split(">")[-1]
except Exception as e:
   print(f"Error extract reactants and product from reaction node: {e}")
</reaction_smiles_extraction>

For traversing the synthetic tree, use this dfs_traverse function:

<schema_dfs_traversal>
def dfs_traverse(node):
    # Process the current node
    # If the node contains children, traverse each recursively
    for child in node.get("children", []):
        dfs_traverse(child)
</schema_dfs_traversal>

Important note: A list of functional group names and reaction type names is appended to the front of the input string. Make sure to account for this in your analysis and modifications.
There is also a list of common RDKIT functions you can use to analyze the code.
This is not a complete set of all RDKIT functions and classes, you may use additional functions and classes as needed.

You have access to the following classes and functions:
- checker.check_fg(fg_name, query_mol_smiles) -> bool
- checker.check_reaction(rxn_type_name, query_rxn_smiles) -> bool
- checker.get_fg_smart(fg_name) -> str
- checker.get_reaction_smart(rxn_type_name) -> str
- checker.get_fg_atom_indices(fg_name, query_mol_smiles) -> list(tuple(tuple(int)))
- checker.check_ring(ring_name, query_mol_smiles) -> bool
- checker.get_ring_smiles(ring_name) -> str
- checker.get_ring_atom_indices(ring_name, query_mol_smiles) -> list(tuple(tuple(int)))
Please perform a rigorous analysis of the code, following these steps:

Follow this Chain of Draft process:

DRAFT 1: Quick Analysis (10-50 words per point)
Think about the following things when conducting your analysis
  - Summarize the intended purpose of the function based on its description.
  - Examine the code line by line, explaining what each part does.
  - List out each variable and function parameter, explaining their purpose.
  - Identify any discrepancies between the code's actual behavior and its intended purpose.
  - Break down the chemical process step-by-step, considering reactants, products, and intermediate steps.
  - Analyze the chemical intuition behind the process, explaining why certain steps or checks are necessary.
  - Identify potential edge cases or unusual chemical scenarios that the code should handle.
  - Remember late-stage corresponds to a low depth count in the DFS traversal, and early-stage corresponds to a high depth count.
  - Remember that the reaction smiles are in the forward reaction direction, but you traverse the tree in the retrosynthetic direction.
  - Ensure that every segment of code accurately reflects the chemical process it is meant to represent. Be wary of any discrepancies.
  - Remember you can use the atom-mapping in the reaction SMILES to track Functional groups and rings between reactants and products.

DRAFT 2: Focused Modifications (5-10 words per change)
For code modifications, keep the following in mind
  - If the code doesn't match the description in a rigerous fashion, suggest specific changes to address these issues.
  - For each proposed change, write out the current code snippet and the proposed modification side by side.
  - Explain the reasoning behind each proposed modification, considering both the technical and chemical aspects.
  - Use checker.check_fg(name, mol_smiles), checker.check_reaction(name, rxn_smiles) and checker.check_ring(name, mol_smiles) functions instead of making up SMARTs patterns.
  - For all reaction types provided to you, always use check_reaction directly. Only try to capture the reaction type by for specific functional groups if the reaction type is not provided.
  - If you are checking for a reaction type included in the set provided to you, assume that if it returns false, the reaction is not present. In this case, do not try to capture the reaction type by checking for specific functional groups.
  - For all functional group and ring patterns provided to you
  - For functional group interconversions, ensure that your code reliably checks that the correct FG is actually changed in the reaction, by concurrently checking for the relevant reaction.
  - Do not modify the function signature; it must always return a single boolean value.
  - Assume all reaction and molecular smiles present in the synthetic route are valid.
- Code change 1: [Before] → [After]
  Reason: [10-50 word explanation]
- Code change 2: [Before] → [After]
  Reason: [10-50 word explanation]
- [Continue as needed]

DRAFT 3: Final Code
Remember to include print statements in the code to aid debugging.
<code>
[Final optimized code]
</code>

Key reminders:
- Low depth = late-stage; high depth = early-stage
- Reaction SMILES = forward direction; traversal = retrosynthetic
- Reaction smiles are atom-mapped where possible.
- Function returns single boolean
- Use checker functions instead of custom SMARTS
- Check correct functional group changes in reactions
- Don't over-engineer, focus on chemical validity
- For sequence checking, ensure correct traversal order
- Function must return a single boolean value

Now analyze the provided code using this Chain of Draft approach.
You must return the code within <code></code> tags. 
"""
