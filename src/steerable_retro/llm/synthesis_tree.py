"""
Evaluate a full route against a query, and demonstrate batch processing of
JSON files in ../data/synth_bench_processed, attaching LLM descriptions.
"""

import asyncio
import base64
import importlib
import json
import os
from io import BytesIO
from typing import Any, Dict, List, Optional

import networkx as nx  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import weave  # type: ignore
from aizynthfinder.chem import FixedRetroReaction  # type: ignore
from aizynthfinder.reactiontree import ReactionTree  # type: ignore
from dotenv import load_dotenv  # type: ignore
from PIL.Image import Image
from prompts import *
from pydantic import BaseModel, model_validator  # type: ignore
from weave.trace.context.call_context import get_current_call  # type: ignore

from steerable_retro.logger import setup_logger
from steerable_retro.utils.rxnimg import get_rxn_img

from .llm_router import router

logger = setup_logger(__name__)


class LM(BaseModel):

    model: str = "gpt-4o"
    vision: bool = False
    prefix: str = ""
    suffix: str = ""
    prompt: Optional[str] = None  # Path to the prompt module
    project_name: str = ""

    async def run(self, tree: ReactionTree, query: str, task: Any):
        """
        Perform an LLM call on the ReactionTree, returning structured data
        about the route's "score" or "description".
        """
        if self.model == "random":
            # For illustration; pick a random score
            response = dict(
                response=f"<score>{np.random.choice(np.arange(1, 11))}</score>",
                url="",
            )
        else:
            rxn_msgs = self.make_msg_sequence(tree)
            response = await self._run_llm(
                rxn_msgs, query, taskid=task.id if task else ""
            )
        return response

    @weave.op()
    async def _run_llm(self, msgs, query, taskid=""):
        """
        Call your LLM endpoint via the router.
        """
        try:
            response = await router.acompletion(
                model=self.model,
                temperature=0.1,
                messages=[
                    {
                        "role": "user",
                        "content": [
                            {
                                "type": "text",
                                "text": self.prefix.format(query=query),
                            },
                            *msgs,
                            {"type": "text", "text": self.suffix},
                        ],
                    },
                ],
            )
            response = response.choices[0].message.content
        except Exception as e:
            logger.error(f"{e}")
            response = "<score>-1</score>"
        current_call = get_current_call()
        return dict(
            response=response,
            url=current_call.ui_url or "-",
        )

    def make_msg_sequence(self, tree: ReactionTree):
        """
        Build a list of "text" or "image_url" messages for each reaction in the route.
        """
        rxns = self.get_smiles_with_depth(tree)
        msgs = []
        for i, (depth, smi) in enumerate(rxns):
            if self.vision:
                inp = self._get_img_msg(smi)
            else:
                inp = self._get_txt_msg(smi)
            msgs.append(
                {"type": "text", "text": f"Reaction #{i+1}. Depth: {depth}"}
            )
            msgs.append(inp)
        return msgs

    def _get_txt_msg(self, smi):
        """Get text message."""
        return {"type": "text", "text": f"{smi}"}

    def _get_img_msg(self, smi):
        """Get encoded reaction image as a data URL."""
        img = get_rxn_img(smi)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        b64img = base64.b64encode(buffered.getvalue()).decode()
        if not b64img:
            raise ValueError("Failed to retrieve the image.")
        msg = {
            "type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{b64img}"},
        }
        return msg

    async def run_single_route(
        self, task, route_dict: Dict, query: str = ""
    ) -> Dict:
        """
        Turn a route dictionary into a ReactionTree, run the LLM, and attach results.
        """
        # task can be None or an object with .id/.prompt
        rt = ReactionTree.from_dict(route_dict)
        result = await self.run(rt, query, task)
        # By default, store the entire LLM response. You could parse it further if needed.
        route_dict["description"] = result[
            "response"
        ]  # The user specifically asked for "description"
        # Also store the raw structured data if you like:
        route_dict["lmdata"] = dict(
            query=query,
            response=result["response"],
            weave_url=result["url"],
            routescore=self._parse_score(result["response"]),
        )
        return route_dict

    async def run_single_task(self, task, data: List[Dict], nroutes=10):
        """
        Demonstration that can run multiple reaction-tree dicts in sequence.
        """
        # limit to nroutes
        data = data[:nroutes]
        results = []
        for d in data:
            updated = await self.run_single_route(
                task, d, task.prompt if task else ""
            )
            results.append(updated)
        return results

    @model_validator(mode="after")
    def load_prompts(self):
        """
        Load environment variables and, if specified, a custom prompt module.
        """
        load_dotenv()
        if self.project_name:
            weave.init(self.project_name)
        if self.prompt is not None:
            module = importlib.import_module(self.prompt)
            self.prefix = module.prefix
            self.suffix = module.suffix
        return self

    @staticmethod
    def _parse_score(response: str) -> float:
        """
        Extract a numeric <score>…</score> from the LLM response, if present.
        """
        try:
            return float(response.split("<score>")[1].split("</score>")[0])
        except Exception:
            return -1  # default if no <score> tag found

    def get_smiles(self, tree: ReactionTree) -> List[str]:
        """Get all reaction SMILES from a tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">>")
                rvsmi = f"{rsmi[1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
        return smiles

    def get_smiles_with_depth(self, tree: ReactionTree) -> List[tuple]:
        """
        Get all reaction SMILES from a tree, along with the depth of each reaction.
        """
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">>")
                rvsmi = f"{rsmi[1]}>>{rsmi[0]}"
                # distance of node m from root
                depth = nx.shortest_path_length(
                    tree.graph, source=tree.root, target=m
                )
                # correct for molecule layers
                depth = int((depth - 1) / 2)
                smiles.append((depth, rvsmi))
        return smiles

    async def process_file(
        self, file_path: str, output_dir: str, query: str = "", top_k: int = 3
    ) -> None:
        """
        Read routes from a JSON file, run LLM on up to top_k routes per target,
        attach 'description', then write the updated file to output_dir.
        """
        logger.info(f"Processing file: {file_path}")
        with open(file_path, "r") as f:
            data = json.load(
                f
            )  # e.g. { "targetA": {"0": {...}, "1": {...}}, "targetB": {...}, ... }

        tasks = []
        # For each target, we randomly sample up to top_k routes
        for target, routes in data.items():
            idxes = list(routes.keys())
            selected_idxs = np.random.choice(
                idxes, min(top_k, len(idxes)), replace=False
            )
            for idx in selected_idxs:
                route_dict = routes[idx]
                # Create an async task that runs run_single_route on route_dict
                tasks.append(
                    asyncio.create_task(
                        self.run_single_route(None, route_dict, query=query)
                    )
                )

        await asyncio.gather(*tasks)

        # After the tasks complete, 'description' is inserted in each route structure.
        os.makedirs(output_dir, exist_ok=True)
        out_file_path = os.path.join(output_dir, os.path.basename(file_path))
        with open(out_file_path, "w") as outf:
            json.dump(data, outf, indent=2)
        logger.info(f"Wrote updated file: {out_file_path}")

    async def process_synth_routes(
        self, input_dir: str, output_dir: str, query: str = "", top_k: int = 3
    ):
        """
        Parallelize file-by-file reading and writing of routes in a directory,
        each generating an LLM-based 'description' field.
        """
        logger.info(f"Reading directory: {input_dir}")

        # Gather all JSON files
        files = [
            f for f in os.listdir(input_dir) if f.lower().endswith(".json")
        ]
        tasks = []
        for fname in files:
            in_path = os.path.join(input_dir, fname)
            tasks.append(
                asyncio.create_task(
                    self.process_file(in_path, output_dir, query, top_k)
                )
            )

        # Run all tasks concurrently
        await asyncio.gather(*tasks)
        logger.info("All files processed.")


async def main():
    """
    Example entrypoint showing how to instantiate LM and run descriptions on a directory.
    """
    # Provide your LLM parameters here
    lm = LM(
        prompt="steer.llm.prompts.fullroute",
        model="claude-3-5-sonnet",
        project_name="steer-test",
    )
    input_path = "../data/synth_bench_processed"
    output_path = "../data/synth_bench_described"
    query_text = "Generate a short description for each route."

    # Process the routes in parallel, file by file
    await lm.process_synth_routes(
        input_path, output_path, query=query_text, top_k=3
    )

    logger.info("Done with route descriptions.")


if __name__ == "__main__":
    asyncio.run(main())
