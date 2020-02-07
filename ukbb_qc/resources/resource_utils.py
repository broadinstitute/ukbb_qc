import argparse
import hail as hl
import logging
from typing import Any, Dict, List, Optional, Tuple, Union
from gnomad_hail.resources.resource_utils import DataException


CURRENT_FREEZE = 5
DATA_SOURCES = ['regeneron', 'broad']
FREEZES = [4, 5]
CURRENT_HAIL_VERSION = "0.2"
