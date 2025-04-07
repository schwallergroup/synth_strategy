"""API for Steerable Retro - create descriptors by running all funcitons."""

import json
from steerable_retro import lm_code
from steerable_retro.logger import setup_logger

logger = setup_logger(__name__)


def run_sequentially(route, *functions):
    """Run all descriptor functions sequentially."""
    results = {}
    for func in functions:
        if callable(func):
            try:
                result = (
                    func(route)
                )  # Call the function, assuming no arguments; adjust if necessary
                results[func.name] = result
                logger.debug(func.description)
            except:
                logger.exception(f"Error running function {func.name}: {func.description}")
                results[func.name] = None
        else:
            print(f"{func} is not a callable function.")
    return results


if __name__ == "__main__":
    # Make a list of all modules in lm_code
    funcs = [
        getattr(lm_code, obj)
        for obj in dir(lm_code)
        if callable(getattr(lm_code, obj)) and obj.startswith("code_")
    ]
    
    # Load test route
    with open("data/test/test_route.json", 'r') as f:
        route = json.load(f)

    results = run_sequentially(route, *funcs)
    logger.info("Results from executed functions:")
    logger.info(results.items())
    for func_name, result in results.items():
        logger.info(f"{func_name}: {result}")
