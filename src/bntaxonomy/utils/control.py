from __future__ import annotations

import contextlib
import json
import os
import re
import sys

from algorecell_types import ReprogrammingStrategies


@contextlib.contextmanager
def suppress_console_output():
    with open(os.devnull, "w") as devnull:
        # suppress stdout and
        orig_stdout_fno = os.dup(sys.stdout.fileno())
        os.dup2(devnull.fileno(), 1)
        orig_stderr_fno = os.dup(sys.stderr.fileno())
        os.dup2(devnull.fileno(), 2)
        try:
            yield
        finally:
            # restore
            os.dup2(orig_stdout_fno, 1)
            os.dup2(orig_stderr_fno, 2)


def check_smaller(p1: dict[str, int], p2: dict[str, int], strict=False):
    is_small = all(p2.get(k, -1) == v for k, v in p1.items())
    if is_small and strict:
        is_small = is_small and not (p1.keys() == p2.keys())
    return is_small


class CtrlResult:
    def __init__(self, name: str, d_list: list[dict[str, int]]) -> None:
        self.name = name
        self.d_list = d_list

    def __repr__(self) -> str:
        return f"CtrlResult({self.name})"

    def __str__(self) -> str:
        return f"{self.d_list}"

    def iter_ctrl_not_included_by(self, other: CtrlResult):
        return (
            x for x in self.d_list if not any(check_smaller(y, x) for y in other.d_list)
        )

    def is_stronger_than(self, other: CtrlResult) -> bool:
        for ctrl in other.iter_ctrl_not_included_by(self):
            return False
        return True

    def dump(self, fname):
        with open(fname, "w") as _f:
            json.dump(self.d_list, _f)

    def sort_d_list(self):
        d_list = [dict(sorted(x.items())) for x in self.d_list]
        d_list.sort(key=lambda x: (len(x), sorted(x.items())))
        self.d_list = d_list

    def drop_nonminimal(self):
        d_list = list()
        for ctrl in self.d_list:
            if not any(True for other in d_list if check_smaller(other, ctrl)):
                d_list.append(ctrl)
        self.d_list = d_list

    def drop_size_limit(self, size_limit: int):
        d_list = list()
        for ctrl in self.d_list:
            if len(ctrl) <= size_limit:
                d_list.append(ctrl)
        self.d_list = d_list


def refine_pert(s: ReprogrammingStrategies):
    text = str(s.perturbations())
    d_list = []
    pattern = re.compile(r"PermanentPerturbation\(([^)]*)\)")
    for match in pattern.findall(text):
        d = dict()
        for item in match.split(", "):
            if "=" not in item:
                continue
            k, v = item.split("=")
            d[k] = int(v)
        d_list.append(d)
    return d_list
