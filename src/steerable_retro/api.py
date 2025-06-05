"""API for Steerable Retro - create descriptors by running all funcitons."""

import json
from typing import List, Callable
from pydantic import BaseModel
from steerable_retro import lm_code
from steerable_retro.logger import setup_logger

logger = setup_logger(__name__)


def run_sequentially(route, *functions):
    """Run all descriptor functions sequentially."""
    results = {}
    for func in functions:
        if callable(func):
            try:
                result = func(route)*1
                results[func.name] = result
                logger.debug(func.description)
            except:
                logger.exception(f"Error running function {func.name}: {func.description}")
                results[func.name] = None
        else:
            print(f"{func} is not a callable function.")
    return results


class StrategyDescriptor(BaseModel):
    functions: List[Callable] = [
        getattr(lm_code, obj)
        for obj in dir(lm_code)
        if callable(getattr(lm_code, obj)) and obj.startswith("code_")
    ]

    def __call__(self, route, *args, **kwargs):
        """Run all descriptor functions sequentially."""
        results = {}
        for func in self.functions:
            if callable(func):
                try:
                    result = func(route, *args, **kwargs)
                    results[func.name] = result * 1.
                    logger.debug(func.description)
                except Exception as e:
                    logger.exception(f"Error running function {func.name}: {func.description}")
                    results[func.name] = None
            else:
                print(f"{func} is not a callable function.")
        return results

    def results_stats(self, results):
        """Calculate some stats for the results."""

        logger.info("Results from executed functions:")
        logger.info(results.items())
        for func_name, result in results.items():
            logger.info(f"{func_name}: {result}")

        # Calculate overall statistics
        total_functions = len(self.functions)
        successful_functions = sum(1 for result in results.values() if result is not None)
        error_functions = total_functions - successful_functions
        positive_balance = sum(i for i in results.values() if i is not None)

        logger.info("Overall statistics:")
        logger.info(f"Total functions: {total_functions}")
        logger.info(f"Successful functions: {successful_functions}")
        logger.info(f"Error functions: {error_functions}")
        logger.info(f"Positive balance: {positive_balance}")




if __name__ == "__main__":
    with open("data/test/test_route.json", 'r') as f:
        route = json.load(f)

    descriptor = StrategyDescriptor()
    results = descriptor(route)
    descriptor.results_stats(results)