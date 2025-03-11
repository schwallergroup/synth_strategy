prefix="""You are an expert computational chemist tasked with analyzing synthetic routes for chemical compounds. Your goal is to identify distinct synthetic strategies that can be implemented using rule-based cheminformatics algorithms. Your analysis will be used to classify thousands of routes into 50-100 distinct strategy groups.

First, carefully review the following synthetic tree:
<synthetic_tree>
"""
suffix="""
</synthetic_tree>

Your task is to analyze this synthetic route with a focus on features that can be programmatically detected. Use the following JSON schema to understand the structure of the synthetic route:

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
rsmi = node["metadata"]["rsmi"]
reactants = rsmi.split(">")[0].split(".")
product = rsmi.split(">")[-1]
</reaction_smiles_extraction>

For traversing the synthetic tree, use this dfs_traverse function:

<schema_dfs_traversal>
def dfs_traverse(node):
    # Process the current node
    # If the node contains children, traverse each recursively
    for child in node.get("children", []):
        dfs_traverse(child)
</schema_dfs_traversal>

In your analysis, focus on these key aspects:

1. Topology changes: Ring formations/openings, chain extensions, scaffold construction approaches
2. Bond disconnection patterns: Key C-C, C-heteroatom bond formation events
3. Fragment combination strategy: Convergent vs. linear, number of fragments merged
4. Functional group handling: Protection/deprotection sequences, functional group interconversions
5. Stereochemistry approach: When and how stereocenters are established
6. Transformation sequence patterns: Order of operations (e.g., reduction before oxidation)

Guidelines for your analysis:
- Identify features that could be detected through substructure matching or reaction classification
- Specify count-based metrics where possible (e.g., "involves 2+ ring formations")
- Note structural motifs that could be searched with SMARTS patterns. For non-aromatic functional groups, provide a name. For aromatic functional groups, provide a SMARTS pattern.
- Focus only on features that cheminformatics tools can reliably detect without requiring chemical intuition

Begin your analysis by breaking down the synthetic route according to these guidelines. Wrap your synthetic route breakdown in <synthetic_route_breakdown> tags. In this breakdown:

1. Identify and list all molecules and reactions in the synthetic tree.
2. For each reaction:
   a. Analyze topology changes (ring formations/openings, chain extensions, etc.)
   b. Identify key bond disconnections (C-C, C-heteroatom)
   c. Note any functional group transformations or protections/deprotections
3. Examine the overall route for:
   a. Fragment combination strategy (convergent vs. linear, number of fragments)
   b. Transformation sequence patterns
4. Note any stereochemistry-related operations
5. Count key features like ring formations, bond disconnections, and fragments
6. Identify specific structural motifs and functional groups that can be detected computationally

Remember to focus only on features that cheminformatics tools can reliably detect without requiring chemical intuition. 
The strategic features you extract should be focused on multi-step sequences where possible, rather than individual reactions.
Ensure the strategies you identify are meaningful, ie functional group presence only matters if it is part of a larger strategy.
After your analysis, provide your output in the following format:

1. Computational Analysis:
<computational_analysis>
- Quantifiable feature 1: [description with specific substructure or transformation]
- Quantifiable feature 2: [description with specific substructure or transformation]
...
</computational_analysis>

2. Extractable Strategy:
<extractable_strategy>
[Concise, rule-implementable description that combines the key features. The strategies should be testable via a binary function that takes a synthetic route as input and returns a boolean value.]
</extractable_strategy>

3. RDKit Considerations:
<rdkit_considerations>
[Notes on how these features could be implemented with RDKIT, SMARTS patterns ( without atom mapping and with [#INT] notation for atoms ), or reaction classification]
</rdkit_considerations>

4. Generated Code:
<code>
[Python functions corresponding to the RDKIT rules extracted from the synthetic routes. These functions should be able to process synthetic routes in the format of the provided JSON schema.
There should be one function for each extracted strategy, which takes the synthetic route, traverses it using the dfs_traverse function, and returns True or False.

Each function should:
1. Include the dfs_traverse call within itself
2. Be self-contained and output a boolean
3. Have a short comment describing which aspect of synthetic strategy it's detecting
4. Include print statements for debugging if needed

Example structure:

def strategy_name(route):
  
    # This function detects [brief description of the strategy].

    def dfs_traverse(node):
        # Strategy-specific logic here
        ...
    
    # Call dfs_traverse on the root node
    dfs_traverse(route)
    
    # Return True or False based on the strategy detection
    return result

]
</code>

Remember that each individual synthetic strategy should be self-contained with a unique function to check it. 
Your goal is to extract patterns that can be implemented in rule-based code to classify thousands of routes into 50-100 distinct strategy groups. Focus only on features that cheminformatics tools can reliably detect without requiring chemical intuition. The functions must simply return True or False if the synthetic strategy is present in the route.
"""