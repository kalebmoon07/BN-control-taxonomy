# stablemotif.py  — cached, registry-friendly runners

import json
import os
import pystablemotifs as sm

from bntaxonomy.iface import register_tool
from bntaxonomy.utils.log import time_check

PRIME_JSON_FILE = "pystablemotif_primes.json"


# -----------------------
# In-memory experiment cache
# -----------------------
# Per expid we keep:
#   - "primes": SM primes dict (from file or JSON cache)
#   - "attrs":  AttractorRepertoire (constructed once per expid)
_cache: dict[int, dict[str, object]] = {}


# -----------------------
# Disk cache helpers
# -----------------------
def _sm_cache_file(cachedir: str) -> str:
    return f"{cachedir}/{PRIME_JSON_FILE}" if cachedir else ""


@time_check
def _try_load_primes(cachedir: str) -> dict | None:
    fname = _sm_cache_file(cachedir)
    if not fname or not os.path.isfile(fname):
        return None
    try:
        with open(fname) as f:
            return json.load(f)
    except Exception:
        return None


@time_check
def _save_primes(primes: dict, cachedir: str) -> None:
    fname = _sm_cache_file(cachedir)
    if not fname:
        return
    try:
        with open(fname, "w") as f:
            json.dump(primes, f)
    except Exception:
        # best-effort; ignore write errors
        pass


# -----------------------
# Preprocessing (cached)
# -----------------------
@time_check
def make_sm_primes_iface(bnet_fname: str, cachedir: str = "") -> dict:
    """
    Load primes from on-disk JSON cache if present; otherwise import from file and cache.
    """
    primes = _try_load_primes(cachedir)
    if primes is None:
        # Import from .bnet (or equivalent) and write cache
        primes = sm.format.import_primes(bnet_fname)
        _save_primes(primes, cachedir)
    return primes


@time_check
def make_sm_attrs_iface(sm_primes: dict):
    """
    Construct an AttractorRepertoire from primes (in-memory only).
    """
    return sm.AttractorRepertoire.from_primes(sm_primes)


# -----------------------
# Runner classes
# -----------------------
class _SM_Base:
    """
    Uniform runner API:
      run(bnet_fname, max_size, target, exclude, expid, cachedir, **kwargs)
      free_experiment(expid)
    - uses_cache: True
    - bn_type: "bnet_file" (expects a path to a Boolean network file)
    """

    uses_cache = True
    bn_type = "bnet_file"

    @classmethod
    def _ensure_primes(cls, expid: int, bnet_fname: str, cachedir: str) -> dict:
        bucket = _cache.setdefault(expid, {})
        if "primes" not in bucket:
            bucket["primes"] = make_sm_primes_iface(bnet_fname, cachedir)
        return bucket["primes"]

    @classmethod
    def _ensure_attrs(cls, expid: int, primes: dict):
        bucket = _cache.setdefault(expid, {})
        if "attrs" not in bucket:
            bucket["attrs"] = make_sm_attrs_iface(primes)
        return bucket["attrs"]

    @staticmethod
    def free_experiment(expid: int):
        _cache.pop(expid, None)


@register_tool
class SM_BruteForce(_SM_Base):
    """
    Driver search via knock_to_partial_state (brute force).
    """

    name = "SM[brute-force]"

    @classmethod
    @time_check
    def run(
        cls,
        bn: str,
        max_size: int,
        target: dict,
        exclude: list,  # kept for API parity
        expid: int,
        cachedir: str,
        **kwargs,
    ):
        primes = cls._ensure_primes(expid, bn, cachedir)
        results = sm.drivers.knock_to_partial_state(
            target, primes, min_drivers=0, max_drivers=max_size, **kwargs
        )
        return results


class _SM_TrapSpaceBase(_SM_Base):
    """
    Trap-space reprogramming via AttractorRepertoire.reprogram_to_trap_spaces.
    - In this package, we only use option `target_method`= `merge` (option `history` gives sequential controls).
    - `driver_method` can be `minimal` or `internal`.

    """

    driver_method: str  # to be set in subclasses

    @classmethod
    @time_check
    def run(
        cls,
        bnet_fname: str,
        max_size: int,
        target: dict,
        exclude: list,  # kept for API parity
        expid: int,
        cachedir: str,
        **kwargs,
    ):
        primes = cls._ensure_primes(expid, bnet_fname, cachedir)
        attrs = cls._ensure_attrs(expid, primes)
        results = attrs.reprogram_to_trap_spaces(
            target,
            target_method="merge",
            driver_method=cls.driver_method,
            max_drivers=max_size,
            **kwargs,
        )
        return results


@register_tool
class SM_TrapSpace_Minimal(_SM_TrapSpaceBase):

    name = "SM[minimal]"
    driver_method = "minimal"


@register_tool
class SM_TrapSpace_Internal(_SM_TrapSpaceBase):

    name = "SM[internal]"
    driver_method = "internal"
