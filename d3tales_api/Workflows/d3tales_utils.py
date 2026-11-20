"""Compatibility import path for the database adapters now under ``D3database``."""

import sys

from d3tales_api.D3database import workflow_adapter

sys.modules[__name__] = workflow_adapter
