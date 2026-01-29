import argparse
import json
import os
from automate_test import insts_ce, insts_practical, get_insts

from bntaxonomy.hierarchy import MultiInputSummary


def main(argv=None):
    cwd = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Summarize counterexamples across instance groups.")
    parser.add_argument(
        "-ig",
        "--inst_groups",
        nargs="+",
        help="Instance-group directories (under 'instances'). If provided, these will be translated to corresponding 'results' folders.",
        default=None,
    )
    parser.add_argument(
        "-i",
        "--instances",
        nargs="+",
        help="Explicit instance folders (under 'instances').",
        default=None,
    )

    args = parser.parse_args(argv)

    # Determine summary source: priority: explicit args > defaults from automate_test
    if args.inst_groups:
        inst_groups = [p.replace("instances", "results") for p in args.inst_groups]
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")
    elif args.instances:
        hc = MultiInputSummary.from_instances(args.instances, "Hierarchy")
    else:
        # default behaviour: use lists from automate_test
        inst_groups = [
            result_dir.replace("instances", "results")
            for result_dir in [insts_practical, insts_ce]
        ]
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")

    # Print counterexample edges and save summary graphs/csvs
    for a1, a2 in hc.ce_G.edges:
        print(
            f"{a1} -> {a2}:",
            [inst.name for inst in hc.ce_G.get_edge_data(a1, a2)["counterexamples"]],
        )

    os.makedirs(f"{cwd}/results", exist_ok=True)
    hc.save(f"{cwd}/results/_summary")

    with open(f"{cwd}/results/counterexample_first_match.csv", "w") as f:
        f.write(hc.to_conflict_matrix_csv(full_ce=False))
    with open(f"{cwd}/results/counterexample_full_match.csv", "w") as f:
        f.write(hc.to_conflict_matrix_csv(full_ce=True))

if __name__ == "__main__":
    main()