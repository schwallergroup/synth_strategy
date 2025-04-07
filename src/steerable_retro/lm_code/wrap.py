from typing import Callable, List, Optional

from pydantic import BaseModel


class LMFunction(BaseModel):
    name: str = ""
    description: str = ""
    func: Callable | None = None

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)
