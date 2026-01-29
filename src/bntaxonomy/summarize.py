import argparse
import json
import os
from bntaxonomy.hierarchy import MultiInputSummary


def main(argv=None):
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
        instances_root = os.path.join("experiments", "instances")
        if not os.path.isdir(instances_root):
            raise FileNotFoundError(f"Instances directory not found: {instances_root}")
        group_dirs = sorted(
            (d for d in os.listdir(instances_root) if os.path.isdir(os.path.join(instances_root, d))), reverse=True
        )
        if not group_dirs:
            raise RuntimeError(f"No instance groups found in {instances_root}")
        inst_groups = [os.path.join("experiments", "results", d) for d in group_dirs]
        hc = MultiInputSummary.from_inst_groups(inst_groups, "Hierarchy")

    # Print counterexample edges and save summary graphs/csvs
    for a1, a2 in hc.ce_G.edges:
        print(
            f"{a1} -> {a2}:",
            [inst.name for inst in hc.ce_G.get_edge_data(a1, a2)["counterexamples"]],
        )

    os.makedirs(f"experiments/results", exist_ok=True)
    hc.save(f"experiments/results/_summary")
    
    exp_names = hc.get_exp_names()
    group_names = hc.get_exp_group_names()
    with open(f"experiments/results/counterexample_group_list.json", "w") as f:
        json.dump(group_names, f, indent=4)
        
    with open(f"experiments/results/counterexample_first_match.csv", "w") as f:
        f.write(hc.to_conflict_matrix_csv(full_ce=False))
    with open(f"experiments/results/counterexample_full_match.csv", "w") as f:
        f.write(hc.to_conflict_matrix_csv(full_ce=True))

if __name__ == "__main__":
    main()