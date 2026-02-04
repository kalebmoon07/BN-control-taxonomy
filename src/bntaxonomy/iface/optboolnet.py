from bntaxonomy.iface import register_tool
from bntaxonomy.utils.log import time_check
from optboolnet.launch import control_fixpoint, control_sync_attr_no_separation


@register_tool
class OptBoolNetFixPoints:
    name = "optbn[FP]"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude):
        s = control_fixpoint(bn, max_size, target, exclude)
        return [item[0] for item in s.perturbations()]


@register_tool
class OptBoolNetSyncAttr:
    name = "optbn[SA]"
    bn_type = "colomoto.BooleanNetwork"

    @time_check
    @staticmethod
    def run(bn, max_size, target, exclude):
        s = control_sync_attr_no_separation(bn, max_size, target, exclude)
        return [item[0] for item in s.perturbations()]
