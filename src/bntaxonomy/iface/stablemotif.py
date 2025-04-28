import pystablemotifs as sm

from bntaxonomy.utils.control import CtrlResult


# preprocessing


def make_sm_primes_iface(bnet_fname: str):
    return sm.format.import_primes(bnet_fname)


def make_sm_attrs_iface(sm_primes: dict):
    return sm.AttractorRepertoire.from_primes(sm_primes)


# control strategies


def ctrl_sm_brute_force_iface(
    sm_primes: dict, target: dict[str, int], max_size: int, **kwargs
):
    results = sm.drivers.knock_to_partial_state(
        target, sm_primes, min_drivers=0, max_drivers=max_size, **kwargs
    )
    return CtrlResult("SM-bf", results)


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
