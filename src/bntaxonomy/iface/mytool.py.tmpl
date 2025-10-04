
"""
Installation instructions

...
"""

from bntaxonomy.iface import register_tool
from bntaxonomy.utils.log import main_logger

cache = {}

def make_primes(bnfile, cachedir):
    from bntaxonomy.iface.pbn import make_pbn_primes_iface, PRIME_JSON_FILE
    if cachedir:
        try:
            with open(f"{cachedir}/{PRIME_JSON_FILE}") as _f:
                primes = json.load(_f)
            main_logger.info("Loaded precomputed pyboolnet primes successfully")
            return primes
        except:
            main_logger.info("Loading precomputed pyboolnet primes fails")
    primes = make_pbn_primes_iface(bnfile)
    if cachedir:
        with open(f"{cachedir}}/{PRIME_JSON_FILE}", "w") as _f:
            json.dump(primes, _f)
    return primes


@register_tool
class MyToolMethod1:
    name = "myGreatMethod"
    uses_cache = True
    bn_type = "bnet_file" # or "colomoto.BooleanNetwork"

    @staticmethod
    def run(bn: str, max_size:int, target:dict, exclude:list,
            expid:int, persistent_cachedir:str):
        if expid not in cache:
            cache[expid] = make_primes(bn, cachedir)
        # must return a list of dict

    @staticmethod
    def free_experiment(expid):
        del cache[expid]
