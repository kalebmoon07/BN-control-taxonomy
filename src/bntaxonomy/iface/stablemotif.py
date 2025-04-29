import pystablemotifs as sm

from bntaxonomy.utils.control import CtrlResult
from bntaxonomy.utils.log import time_check

PRIME_JSON_FILE = "pystablemotif_primes.json"


# preprocessing


@time_check
@DeprecationWarning  # use pyboolnet instead
def make_sm_primes_iface(bnet_fname: str):
    print("Generating primes from pystablemotifs")
    return sm.format.import_primes(bnet_fname)


@time_check
def make_sm_attrs_iface(sm_primes: dict):
    print("Generating AttractorRepertoire from pystablemotifs")
    return sm.AttractorRepertoire.from_primes(sm_primes)


# control strategies


@time_check
def ctrl_sm_brute_force_iface(
    sm_primes: dict, target: dict[str, int], max_size: int, **kwargs
):
    results = sm.drivers.knock_to_partial_state(
        target, sm_primes, min_drivers=0, max_drivers=max_size, **kwargs
    )
    return CtrlResult("SM-bf", results)


@time_check
def ctrl_sm_trapspace_iface(
    sm_attrs: sm.AttractorRepertoire,
    target: dict[str, int],
    target_method: str,
    driver_method: str,
    max_size: int,
    **kwargs,
):
    results = sm_attrs.reprogram_to_trap_spaces(
        target, target_method, driver_method, max_drivers=max_size, **kwargs
    )
    return CtrlResult(f"SM-{target_method}-{driver_method[:3]}", results)
