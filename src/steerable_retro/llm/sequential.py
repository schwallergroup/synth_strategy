"""Evaluate a full route against a query."""

import asyncio
import aiofiles
import base64
import importlib
import os
import json
from io import BytesIO
from typing import Any, Dict, List, Optional


from concurrent.futures import ThreadPoolExecutor, as_completed
import networkx as nx  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import weave  # type: ignore
from aizynthfinder.chem import FixedRetroReaction, RetroReaction  # type: ignore
from aizynthfinder.reactiontree import ReactionTree  # type: ignore
from dotenv import load_dotenv  # type: ignore
from PIL.Image import Image
from pydantic import BaseModel, model_validator  # type: ignore
from weave.trace.context.call_context import get_current_call  # type: ignore

from steer.logger import setup_logger
from steer.utils.rxnimg import get_rxn_img

from .llm_router import router
from .prompts import *

logger = setup_logger(__name__)

weave.init('liac/starting_material_prompt')
class LM(BaseModel):
    """LLM Heuristic for scoring reactions."""

    model: str = "deepseek-r1"
    vision: bool = False
    prefix: str = ""
    suffix: str = ""
    prompt: Optional[str] = None  # Path to the prompt module
    project_name: str = ""

    async def run(self, tree: ReactionTree, query: str):
        """Get smiles and run LLM."""

        if self.model == "random":
            response = dict(
                response=f"<score>{np.random.choice(np.arange(1,11))}</score>",
                url="",
            )
        else:
            rxn_msgs = self.make_msg_sequence(tree)
            response = await self._run_llm(rxn_msgs, query)
        return response

    @weave.op()
    async def _run_llm(self, msgs, query, taskid=""):
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
                                "text": self.prefix
                            },
                            *msgs,
                            {"type": "text", "text": self.suffix},
                        ],
                    },
                ],
            )
            response = response.choices[0].message.content
        except Exception as e:
            print(f"Error running LLM: {e}")
            logger.error(f"{e}")
            response = "<score>-1</score>"

        current_call = get_current_call()
        return dict(
            response=response,
            url=current_call.ui_url or "-",
        )

    def make_msg_sequence(self, tree: ReactionTree):
        rxns = self.get_smiles_with_depth(tree)

        msgs = []
        for i, s in enumerate(rxns):
            depth, smi = s
            if self.vision:
                inp = self._get_img_msg(smi)
            else:
                inp = self._get_txt_msg(smi)

            msg = [
                {"type": "text", "text": f"Reaction #{i+1}. Depth: {depth}"},
                inp,
            ]
            msgs.extend(msg)
        return msgs

    def _get_txt_msg(self, smi):
        """Get text message."""
        return {"type": "text", "text": f"{smi}"}

    def _get_img_msg(self, smi):
        """Get image message."""

        img = get_rxn_img(smi)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        b64img = base64.b64encode(buffered.getvalue()).decode()
        if b64img is None:
            raise ValueError("Failed to retrieve the image.")

        msg = {
            "type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{b64img}"},
        }
        return msg

    async def run_single_route(self, task, d):
        result = await self.run(ReactionTree.from_dict(d), task.prompt, task)
        d["lmdata"] = dict(
            query=task.prompt,
            response=result["response"],
            weave_url=result["url"],
            routescore=self._parse_score(result["response"]),
        )
        return d

    async def run_single_task(self, task, data, nroutes=10):
        result = await asyncio.gather(
            *[self.run_single_route(task, d) for d in data[:nroutes]]
        )
        return result

    @model_validator(mode="after")
    def load_prompts(self):
        load_dotenv()
        if self.project_name:
            weave.init(self.project_name)
        if self.prompt is not None:
            module = importlib.import_module(self.prompt)
            self.prefix = module.prefix
            self.suffix = module.suffix
        return self

    @staticmethod
    def _parse_score(response):
        try:
            return float(response.split("<score>")[1].split("</score>")[0])
        except:
            return -1  # Default score (min)

    def get_smiles(self, tree: ReactionTree):
        """Get all smiles from a tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
        return smiles

    def get_smiles_with_depth(self, tree: ReactionTree):
        """Get all smiles from a tree, with depth in tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"

                # Get distance of node m from root
                depth = nx.shortest_path_length(
                    tree.graph, source=tree.root, target=m
                )
                depth = int(
                    (depth - 1) / 2
                )  # Correct for molecule nodes in between
                smiles.append((depth, rvsmi))
        return smiles


async def process_file(file_path, lm):
    results = []
    hash = file_path.split("_")[-1].replace(".json", "")

    query = ""
                                  
    async with aiofiles.open(file_path, "r") as f:
        contents = await f.read()
        data = json.loads(contents)
    for d in data:
        result = await lm.run(ReactionTree.from_dict(d), f"{query}")
        results.append(result)
    return results

async def main():
    model_aliases = [
        "deepseek-v3",
        # "gpt-4o-2024-11-20",
        # "gpt-4o-2024-05-13",
        # "gpt-4-0314", 
        # "gpt-4-1106-preview",
        # "gpt-4-0613-preview",
        # "gpt-4-turbo-2024-04-09",
    ]
    for model_name in model_aliases:
        lm = LM(
            prompt="steer.llm.prompts.route_desc_sm_fgs",
            model=model_name,
            project_name="route_desc_sm_fgs",
        )
        synth_bench_dir = "../data/starting_material/fg_test"
        output_dir = "../data/starting_material/route_fg"
        print(f"Processing output in {output_dir}")
        # Get list of JSON files in the synth_bench directory
        json_files = [
            os.path.join(synth_bench_dir, f)
            for f in os.listdir(synth_bench_dir)
            if f.endswith('ca06156bee8f14dcf0bd7e14f68eddcc.json')
        ]
        # json_files = ["/home/dparm/reaction_utils/rxnutils/data/pa_routes/ref_routes_n1.json"]
        # output_dir = "/home/dparm/reaction_utils/rxnutils/data/pa_routes"
        # Limit the number of concurrent tasks to prevent resource exhaustion
        semaphore = asyncio.Semaphore(20)  # Adjust the limit as needed

        async def sem_process(file_path):
            async with semaphore:
                try:
                    results = await process_file(file_path, lm)
                    output_file = os.path.join(
                        output_dir,
                        os.path.basename(file_path).replace(
                            ".json", f"{model_name}.json"
                        ),
                    )
                    # Use asynchronous file I/O
                    async with aiofiles.open(output_file, "w") as f:
                        await f.write(json.dumps(results))
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")

        # Create a list of tasks
        tasks = [asyncio.create_task(sem_process(fp)) for fp in json_files]
        
        # Await all tasks
        await asyncio.gather(*tasks)



if __name__ == "__main__":
    asyncio.run(main())
