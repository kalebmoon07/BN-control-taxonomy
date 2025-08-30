from bntaxonomy.iface import register_tool
from bntaxonomy.utils.log import main_logger

from bntaxonomy.utils.control import suppress_console_output
from bntaxonomy.utils.log import time_check

import bonesis
from colomoto.minibn import BooleanNetwork
from bonesis.reprogramming import (
    marker_reprogramming_fixpoints,
    marker_reprogramming,
)

@register_tool
class BoNesisFixedPoints:
    name = "BoNesis[FP]"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude):
        with suppress_console_output():
            return list(marker_reprogramming_fixpoints(bn, target, max_size,
                                       exclude=exclude, at_least_one=False))

@register_tool
class BoNesisTrapSpaces:
    name = "BoNesis[MTS]"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude):
        with suppress_console_output():
            return list(marker_reprogramming(bn, target, max_size,
                                             exclude=exclude))
