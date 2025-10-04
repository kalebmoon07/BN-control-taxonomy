import json
import logging
import os

from bntaxonomy.iface import register_tool
from bntaxonomy.utils.log import time_check, main_logger

from pyboolnet.file_exchange import bnet2primes
from itertools import combinations, product
from typing import List, Optional

from pyboolnet.prime_implicants import find_inputs, find_constants, create_constants, percolate, remove_variables
from pyboolnet.trap_spaces import compute_trap_spaces
from pyboolnet.attractors import completeness
from pyboolnet.model_checking import model_checking
from pyboolnet.temporal_logic import subspace2proposition
from pyboolnet.helpers import dicts_are_consistent

log = logging.getLogger(__file__)


def is_included_in_subspace(subspace1: dict, subspace2: dict) -> bool:
    """
    Test whether *subspace1* is contained in *subspace2*.

    **arguments**:
        * *subspace1*, *subspace2* (dicts): subspaces.
    **returns**:
        * Answer (bool): whether *subspace1* is contained in *subspace2*.
    **example**::
        >>> is_included_in_subspace({'v1': 0, 'v2': 1}, {'v2': 1})
        True
    """

    answer = all(x in subspace1 and subspace1[x] == subspace2[x] for x in subspace2.keys())

    return answer


def EFAG_set_of_subspaces(primes: dict, subspaces: List[dict]) -> str:
    """
    Construct a CTL formula that queries whether there is a path that leads to the union of the *subspaces* and stays there.

    **arguments**:
        * *primes*: prime implicants.
        * *subspaces*: list of subspaces.

    **returns**:
        * *formula*: the CTL formula.
    """

    formula = 'EF(AG(' + ' | '.join([subspace2proposition(primes, x) for x in subspaces]) + '))'

    return formula


def fix_components_and_reduce(primes: dict, subspace: dict, keep_vars: List[str] = []) -> dict:
    """
    Fix the variables fixed in *subspace* and percolates their values in *primes*. Returns the resulting set of primes after removing all the constant variables that are not in *keep_vars*.

    **arguments**:
        * *primes*: prime implicants.
        * *subspace*: subspace.
        * *keep_vars*: list of variables to keep in *primes*.

    **returns**:
        * *new_primes*: prime implicants after fixing the variables in *subspace*, percolating them and removing the constant variables not in *keep_vars*.
    """

    new_primes = percolate(primes, add_constants=subspace, copy=True)
    removable_vars = [k for k in find_constants(new_primes) if k not in keep_vars]
    new_primes = remove_variables(new_primes, removable_vars, copy=True)

    return new_primes


def select_control_strategies_by_percolation(primes: dict, strategies: List[dict], target: List[dict]) -> List[dict]:
    """
    Select the elements of *strategies* that are control strategies for *target* by direct percolation.

    **arguments**:
        * *primes*: prime implicants.
        * *strategies*: list of control strategies.
        * *target*: list of subspaces defining the target subset.

    **returns**:
        * *selected_strategies*: list of control strategies by direct percolation.
    """

    selected_strategies = []
    for x in strategies:
        percolation = find_constants(primes=percolate(primes=primes, add_constants=x, copy=True))
        if any(is_included_in_subspace(percolation, subs) for subs in target):
            selected_strategies.append(x)

    return selected_strategies


def control_is_valid_in_trap_spaces(primes: dict, trap_spaces: List[dict], target: List[dict], update: str) -> bool:
    """
    Return whether the *trap_spaces* are compatible with the *target* by checking that no trap space is disjoint from the *target* and applying the control query to all the trap spaces oscillating in and out the *target*.

    **arguments**:
        * *primes*: prime implicants.
        * *trap_spaces*: list of trap spaces.
        * *target*: list of subspaces defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.

    **returns**:
        * *cs_perc*: list of control strategies by direct percolation.
    """

    if not all(any(dicts_are_consistent(ts, subs) for subs in target) for ts in trap_spaces):
        return False

    half_ts = [ts for ts in trap_spaces if not any(is_included_in_subspace(ts, subs) for subs in target)]
    for ts in half_ts:
        if not reduce_and_run_control_query(primes, ts, target, update):
            return False

    return True


def reduce_and_run_control_query(primes: dict, subspace: dict, target: List[dict], update: str):
    """
    Run the model checking query for control after reducing the network by percolating the values in *subspace*.

    **arguments**:
        * *primes*: prime implicants.
        * *subspace*: subspace.
        * *target*: list of subspaces defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.

    **returns**:
        * True if the query is true, false otherwise.
    """

    target_vars = list(set([item for subs in target for item in subs]))
    new_primes = fix_components_and_reduce(primes, subspace, keep_vars=target_vars)
    answer = run_control_query(new_primes, target, update)

    return answer


def run_control_query(primes: dict, target: List[dict], update: str) -> bool:
    """
    Run the model checking query for control as described in :ref:`CifuentesFontanals2022 <CifuentesFontanals2022>` Sec 4.2.

    **arguments**:
        * *primes*: prime implicants.
        * *target*: list of subspaces defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.

    **returns**:
        * True if the query is true, False otherwise.
    """

    spec = 'CTLSPEC ' + EFAG_set_of_subspaces(primes, target)
    init = "INIT TRUE"
    answer = model_checking(primes, update, init, spec)

    return answer


def control_direct_percolation(primes: dict, candidate: dict, target: List[dict]) -> bool:
    """
    Check whether the subspace *candidate* is a control strategy for *target* by direct percolation.

    **arguments**:
        * *primes*: prime implicants.
        * *candidate*: subspace.
        * *target*: list of subspaces defining the target subset.

    **returns**:
        * True if the *candidate* percolates into the *target*, False otherwise.
    """

    perc = find_constants(primes=percolate(primes=primes, add_constants=candidate, copy=True))

    if any(is_included_in_subspace(perc, subs) for subs in target):
        log.info(f"Intervention (only percolation): {candidate}")
        return True

    return False


def control_completeness(primes: dict, candidate: dict, target: dict, update: str) -> Optional[bool]:
    """
    Check whether the subspace *candidate* is a control strategy for *target* by the completeness approach,
    described in :ref:`CifuentesFontanals2022 <CifuentesFontanals2022>` Sec 3.2.

    **arguments**:
        * *primes*: prime implicants.
        * *candidate*: subspace.
        * *target*: subspace defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.

    **returns**:
        * True if the *candidate* is a control strategy by completeness for *target*, False otherwise.
    """

    if type(target) != dict:
        log.error("The target must be a subspace (dict).")
        return

    perc = find_constants(primes=percolate(primes=primes, add_constants=candidate, copy=True))
    new_primes = fix_components_and_reduce(primes, perc, keep_vars=list(target.keys()))
    minimal_trap_spaces = compute_trap_spaces(new_primes, "min")

    if not all(is_included_in_subspace(T, target) for T in minimal_trap_spaces):
        return False

    if completeness(new_primes, update):
        log.info(f"Intervention (by completeness): {candidate}")
        return True

    return False


def control_model_checking(primes: dict, candidate: dict, target: List[dict], update: str, max_output_trapspaces: int = 10000000) -> Optional[bool]:
    """
    Check whether the subspace *candidate* is a control strategy for *target* using the model checking approach
    described in :ref:`CifuentesFontanals2022 <CifuentesFontanals2022>` Sec 4.3.

    **arguments**:
        * *primes*: prime implicants.
        * *candidate*: subspace.
        * *target*: list of subspaces defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.

    **returns**:
        * True if the *candidate* is a control strategy by completeness for *target*, False otherwise.
    """

    if type(target) != list:
        log.error("The target must be a list of subspaces.")
        return

    perc = find_constants(primes=percolate(primes=primes, add_constants=candidate, copy=True))
    new_primes = fix_components_and_reduce(primes, perc, keep_vars=list(set([item for subs in target for item in subs])))
    minimal_trap_spaces = compute_trap_spaces(new_primes, "min", max_output=max_output_trapspaces)

    if not control_is_valid_in_trap_spaces(new_primes, minimal_trap_spaces, target, update):
        return False

    if run_control_query(new_primes, target, update):
        log.info(f"Intervention (by CTL formula): {candidate}")
        return True

    return False


def find_necessary_interventions(primes: dict, target: List[dict]) -> dict:
    """
    Find the names and values of the inputs and constants from *primes* that are fixed in the *target*.

    **arguments**:
        * *primes*: prime implicants
        * *target*: list of subspaces

    **returns**:
        * *selected_vars*: Names and values that are inputs or constants in *primes* and are fixed in the *target*.
    """

    selected_vars = dict()
    candidates = find_inputs(primes) + list(find_constants(primes).keys())
    for x in candidates:
        if all(x in y.keys() for y in target):
            if all(y[x] == z[x] for y in target for z in target):
                selected_vars[x] = target[0][x]

    return selected_vars


def find_common_variables_in_control_strategies(primes: dict, target: List[dict]) -> dict:
    """
    Find the names and values of the constants from *primes* that are fixed to a different value in the *target* and therefore need to be part of any control strategy.

    **arguments**:
        * *primes*: prime implicants
        * *target*: list of subspaces

    **returns**:
        * Names and values that are constants in *primes* and are fixed to a different value in the *target*.
    """

    common_inputs_and_constants_in_target = find_necessary_interventions(primes, target)
    constants = find_constants(primes)
    right_constants = [x for x in constants if x in common_inputs_and_constants_in_target.keys() and common_inputs_and_constants_in_target[x] == constants[x]]
    common_variables = {k: common_inputs_and_constants_in_target[k] for k in common_inputs_and_constants_in_target.keys() if k not in right_constants}

    return common_variables


def is_control_strategy(primes: dict, candidate: dict, target: List[dict], update: str, max_output: int = 1000000) -> Optional[bool]:
    """
    Check whether the *candidate* subspace is a control strategy for the *target* subset,
    as defined in ref:`CifuentesFontanals2022 <CifuentesFontanals2022>`.

    **arguments**:
        * *primes*: prime implicants.
        * *candidate*: candidate subspace.
        * *target*: list of subspaces defining the target subset.
        * *update*: type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.
        * *max_output*: maximum number of trap spaces computed. Default value: 1000000.

    **returns**:
        * *answer*: True if the *candidate* is a control strategy for the *target*, False otherwise. When False, a counterexample (minimal trap space disjoint with target or state not satisfying the CTL query) is also returned.

    **example**::
        >>> is_control_strategy(primes, {'v1': 1}, [{'v2': 0}, {'v3':1}], "asynchronous")
    """

    if type(target) != list:
        log.error("The target must be a list of subspaces.")
        return

    if control_direct_percolation(primes, candidate, target):
        return True

    new_primes = fix_components_and_reduce(primes, candidate, list(set([x for subs in target for x in subs])))
    minimal_trap_spaces = compute_trap_spaces(new_primes, "min", max_output=max_output)

    if not control_is_valid_in_trap_spaces(new_primes, minimal_trap_spaces, target, update):
        return False

    answer = run_control_query(new_primes, target, update)

    return answer


def compute_control_strategies_with_completeness(primes: dict, target: dict, update: str = "asynchronous", limit: int = 3, avoid_nodes: List[str] = None, starting_length: int = 0, known_strategies: List[dict] = None) -> Optional[List[dict]]:
    """
    Identify control strategies for the *target* subspace using the completeness approach
    described in :ref:`CifuentesFontanals2022 <CifuentesFontanals2022>` Sec 3.2.

    **arguments**:
        *primes*: prime implicants.
        *target*: target subspace.
        *update*: the type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*.
        *limit*: maximal size of the control strategies. Default value: 3.
        *starting_length*: minimum possible size of the control strategies. Default value: 0.
        *known_strategies*: list of already identified control strategies. Default value: empty list.
        *avoid_nodes*: list of nodes that cannot be part of any control strategy. Default value: empty list.

    **returns**:
        * *list_strategies*: list of control strategies for the *target* subspace obtained using completeness.

    **example**::
        >>> control_strategies_completeness(primes, {'v1': 1}, "asynchronous")
    """

    if type(target) == list:
        log.error("The target must be a subspace.")
        return

    avoid_nodes = avoid_nodes or []
    known_strategies = known_strategies or []

    candidate_variables = [x for x in primes.keys() if x not in avoid_nodes]
    list_strategies = known_strategies
    perc_true = known_strategies
    perc_false = []

    common_vars_in_cs = find_common_variables_in_control_strategies(primes, [target])
    candidate_variables = [x for x in primes.keys() if x not in common_vars_in_cs.keys() and x not in avoid_nodes]
    log.info(f"Number of common variables in the CS: {len(common_vars_in_cs)}")
    log.info(f"Number of candiadate variables: {len(candidate_variables)}")


    for i in range(max(0, starting_length - len(common_vars_in_cs)), limit + 1 - len(common_vars_in_cs)):

        log.info(f"Checking control strategies of size {i + len(common_vars_in_cs)}")

        for vs in combinations(candidate_variables, i):

            subsets = product(*[(0, 1)]*i)
            for ss in subsets:
                candidate = dict(zip(vs, ss))
                candidate.update(common_vars_in_cs)

                if not any(is_included_in_subspace(candidate, x) for x in list_strategies):
                    perc = find_constants(primes=percolate(primes=primes, add_constants=candidate, copy=True))

                    if perc in perc_true:
                        log.info(f"Intervention: {candidate}")
                        list_strategies.append(candidate)

                    elif perc not in perc_false:

                        if control_direct_percolation(primes, candidate, [target]):
                            perc_true.append(perc)
                            list_strategies.append(candidate)

                        elif control_completeness(primes, candidate, target, update):
                            perc_true.append(perc)
                            list_strategies.append(candidate)

                        else:
                            perc_false.append(perc)

    return list_strategies


def compute_control_strategies_with_model_checking(primes: dict, target: List[dict], update: str = "asynchronous", limit: int = 3, avoid_nodes: List[str] = None, max_output_trapspaces: int = 1000000, starting_length: int = 0, known_strategies: List[dict] = None) -> Optional[List[dict]]:
    """
    Identify all minimal control strategies for the *target* subset using the model checking approach
    described in :ref:`CifuentesFontanals2022 <CifuentesFontanals2022>` Sec 4.3.

    **arguments**:
        *primes*: prime implicants
        *target*: list of subspaces defining the target subset
        *update*: the type of update, either *"synchronous"*, *"asynchronous"* or *"mixed"*
        *limit*: maximal size of the control strategies. Default value: 3.
        *starting_length*: minimum possible size of the control strategies. Default value: 0.
        *known_strategies*: list of already identified control strategies. Default value: empty list.
        *avoid_nodes*: list of nodes that cannot be part of the control strategies. Default value: empty list.

    **returns**:
        * *list_strategies*: list of control strategies (dict) of *subspace* obtained using completeness.

    **example**::
        >>> ultimate_control_multispace(primes, {'v1': 1}, "asynchronous")

    """

    if type(target) != list:
        log.error("The target must be a list.")
        return

    avoid_nodes = avoid_nodes or []
    known_strategies = known_strategies or []

    list_strategies = known_strategies
    perc_true = known_strategies
    perc_false = []

    common_vars_in_cs = find_common_variables_in_control_strategies(primes, target)
    candidate_variables = [x for x in primes.keys() if x not in common_vars_in_cs.keys() and x not in avoid_nodes]
    log.info(f"Number of common variables in the CS: {len(common_vars_in_cs)}")
    log.info(f"Number of candiadate variables: {len(candidate_variables)}")


    for i in range(max(0, starting_length - len(common_vars_in_cs)), limit + 1 - len(common_vars_in_cs)):

        log.info(f"Checking control strategies of size {i + len(common_vars_in_cs)}")

        for vs in combinations(candidate_variables, i):
            subsets = product(*[(0, 1)]*i)

            for ss in subsets:
                candidate = dict(zip(vs, ss))
                candidate.update(common_vars_in_cs)

                if not any(is_included_in_subspace(candidate, x) for x in list_strategies):
                    perc = find_constants(primes=percolate(primes=primes, add_constants=candidate, copy=True))

                    if perc in perc_true:
                        log.info(f"Intervention: {candidate}")
                        list_strategies.append(candidate)

                    elif perc not in perc_false:

                        if control_direct_percolation(primes, candidate, target):
                            perc_true.append(perc)
                            list_strategies.append(candidate)

                        elif control_model_checking(primes, candidate, target, update):
                            perc_true.append(perc)
                            list_strategies.append(candidate)

                        else:
                            perc_false.append(perc)

    return list_strategies



from clingo import Control
from pyboolnet.trap_spaces import compute_trapspaces_that_intersect_subspace

def run_node_edge_control_asp(program_instance: str):
    
    ctl = Control(arguments=[f"--models=0", "--opt-mode=optN", "--enum-mode=domRec", "--heuristic=Domain", "--dom-mod=5,16",])
    
    ctl.add(name="base", parameters={}, program=program_instance)
    
    ctl.add(name="base", parameters={}, program="""
        goal(T,S) :- goal(Z,T,S), Z < 0.
        satisfy(V,W,S) :- formula(W,D); dnf(D,C); clause(C,V,S).
        closure(V,T)   :- goal(V,T).
        closure(V,S*T) :- closure(W,T); satisfy(V,W,S); not goal(V,-S*T).
        { node(V,S) } :- closure(V,S), not avoid_node(V), satisfied(Z), Z < 0.
        { node(V,S) : goal(Z,V,S), not avoid_node(V), satisfied(Z), subspace(Z), Z >= 0}.
        { edge(Vi,Vj,1); edge(Vi,Vj,-1) } :- formula(Vj,D), dnf(D,C), clause(C, Vi, S), not avoid_edge(Vi,Vj), satisfied(Z), Z < 0.
        { edge(Vi,Vj,1); edge(Vi,Vj,-1) } :- formula(Vj,D), dnf(D,C), clause(C, Vi, S), not avoid_edge(Vi,Vj), goal(Z,Vj,T), satisfied(Z), subspace(Z), Z >= 0.
        :- node(V,S), node(V,-S).
        :- edge(Vi,Vj,S), edge(Vi, Vj, -S).
        :- node(V,S), edge(V,Vj).
        :- node(V), edge(Vi,V).
        node(V) :- node(V,S).
        edge(Vi,Vj) :- edge(Vi,Vj,S).
        new_clause(C,V,S) :- clause(C,V,S); dnf(D,C); formula(Vj,D); not edge(V,Vj).
        remove_dnf(D,C):- clause(C,Vi,-S); edge(Vi,Vj,S); dnf(D,C); formula(Vj,D).
        new_dnf(D,C) :- new_clause(C,Vi,S); dnf(D,C); formula(Vj,D); not remove_dnf(D,C).
        remove_formula(Vj,D) :- dnf(D,C); formula(Vj,D); edge(Vi,Vj,S) : clause(C,Vi,S).
        new_formula(V,D) :- new_dnf(D,C); formula(V,D); not remove_formula(V,D).
        fixed_node(V,1) :- remove_formula(V,D).
        fixed_node(V,-1) :- not remove_formula(V,D); not new_formula(V,D); formula(V,D).
        intervention(V,S) :- node(V,S).
        intervention(V,S) :- not node(V,S), not node(V,-S), fixed_node(V,S).
        intervention(V) :- intervention(V,S).
        eval_formula(Z,V,S) :- subspace(Z); intervention(V,S).
        free(Z,V,D) :- new_formula(V,D); subspace(Z); not intervention(V).
        eval_clause(Z,C,-1) :- new_clause(C,V,S); eval_formula(Z,V,-S).
        eval_formula(Z,V, 1) :- free(Z,V,D); eval_formula(Z,W,T) : new_clause(C,W,T); new_dnf(D,C).
        eval_formula(Z,V,-1) :- free(Z,V,D); eval_clause(Z,C,-1) : new_dnf(D,C).
        not satisfied(Z) :- goal(Z,T,S), not eval_formula(Z,T,S), subspace(Z).
        satisfied(Z) :- eval_formula(Z,T,S) : goal(Z,T,S); subspace(Z).
        0 < { satisfied(Z) : subspace(Z) }.
        :- maxsize>0; maxsize + 1 { node(V,R); edge(Vi,Vj,S) }.
        :- maxnodes<0; 1 { node(V,S) }.
        :- maxedges<0; 1 { edge(Vi,Vj,S) }.
        #show node/2.
        #show edge/3.
        """)
    
    ctl.ground([("base", [])])
    
    models = []
    with ctl.solve(yield_=True) as handle:
        for model in handle:
            models.append(model.symbols(shown=True))

    return models


def read_asp_output(primes: dict, models: List[list]):

    lower_to_prime = {n.lower(): n for n in primes}
    value_to_boolean = {1:1, -1:0}
    
    cs_total = []
    for x in models:
        cs = {}
        for y in x:
            if y.name == "node":
                cs[lower_to_prime[y.arguments[0].name]] = value_to_boolean[y.arguments[1].number]
            if y.name == "edge":
                cs[(lower_to_prime[y.arguments[0].name], lower_to_prime[y.arguments[1].name])] = value_to_boolean[y.arguments[2].number]
        cs_total.append(cs)
    
    return cs_total


def run_control_problem(primes, target, intervention_type, control_type, avoid_nodes: dict = {}, avoid_edges: dict = {}, limit: int = 3, output_file: str = "", use_attractors: bool = True, complex_attractors: List[List[dict]] = []):

    # Setting targets and computing selected trap spaces

    if control_type in ["trap_spaces", "transient", "both"]:
        tspaces = compute_trapspaces_that_intersect_subspace(primes=primes, subspace=target, type_="percolated", max_output=1000000)
        if {} not in tspaces:
            tspaces.append({})
        tsmin = compute_trap_spaces(primes, "min")
        target_trap_spaces = select_trapspaces(tspaces=tspaces, subspace=target, use_attractors=use_attractors, tsmin=tsmin, complex_attractors=complex_attractors)
        target_percolation = []
        print("Num of selected trap spaces:", len(target_trap_spaces))

        if control_type == "transient":
            target_percolation = target_trap_spaces
            target_trap_spaces = []
        elif control_type == "both":
            target_percolation = [target]
    else:
        target_trap_spaces = []
        target_percolation = [target]

    # Computing CS in ASP

    program_instance = create_asp_program_instance(primes=primes, intervention_type=intervention_type, target_trap_spaces=target_trap_spaces, target_subspaces=target_percolation, max_size=limit, avoid_nodes=avoid_nodes, avoid_edges=avoid_edges, filename="program_instance")
    models = run_node_edge_control_asp(program_instance)
    cs_asp = read_asp_output(primes, models)

    # Saving output

    if output_file != "":
        with open(output_file+".py", "w") as f:
            f.write("Target = " + str(target) + "\n")
            f.write("#Control strategies using " + str(control_type) + "\n cs = " + str(cs_asp) + "\n")

    return cs_asp


def create_asp_program_instance(primes: dict, intervention_type: str, target_trap_spaces: List[dict] = [], target_subspaces: List[dict] = [], max_size: int = 3, avoid_nodes: List[str] = [], avoid_edges: List[str] = [], filename: str = "") -> str:
    """
    Encodes the control strategy problem is ASP.
    The output is a string. If *filename* is provided, it saves the output on a file.
    """

    # Encoding of the Boolean function
    nodes_to_avoid = ""
    edges_to_avoid = ""
    for x in avoid_edges:
        edges_to_avoid = edges_to_avoid + "avoid_edge(" + x[0] + "," + x[1] + "). "
    formulas = ""
    dnfs = ""
    clauses = ""
    clauses_dict = dict()
    primes_included = dict()
    id_form = -1
    cont_clause = 0
    for x in primes.keys():
        id_form = id_form + 1
        if x in avoid_nodes:
            nodes_to_avoid = nodes_to_avoid + "avoid_node(" + x + "). "
        formulas = formulas + "formula(" + x + ", " + str(id_form) + "). "
        for p in primes[x][1]:
            id_clause = str(cont_clause)
            cont_clause = cont_clause + 1
            text_clause = ""
            for y in p.keys():
                value = str(p[y])
                if str(p[y]) == "0":
                    value = "-1"
                text_clause = text_clause + "clause(" + id_clause + ", " + y + ", " + value + "). "
            clauses_dict[id_clause] = text_clause
            clauses = clauses + text_clause
            dnfs = dnfs + "dnf(" + str(id_form) + ", " + id_clause + "). "
    # Encoding goal subspaces
    subspaces = ""
    goals = ""

    id_subspace = 0
    for s in target_subspaces:
        id_subspace = id_subspace - 1
        subspaces = subspaces + f"subspace({id_subspace}). "
        for x in s.keys():
            value = str(s[x])
            if str(s[x]) == "0":
                value = "-1"
            goals = goals + "goal(" + str(id_subspace) + ", " + x + ", " + value + "). "

    id_subspace = -1
    for s in target_trap_spaces:
        id_subspace = id_subspace + 1
        subspaces = subspaces + f"subspace({id_subspace}). "
        for x in s.keys():
            value = "-1" if str(s[x]) == "0" else str(s[x])
            goals = goals + "goal(" + str(id_subspace) + ", " + x + ", " + value + "). "

    max_nodes = max_size
    max_edges = max_size
    if intervention_type == "node":
        max_edges = -1
    if intervention_type == "edge":
        max_nodes = -1
    constants = "#const maxsize=" + str(max_size) + "." + "\n\n" + "#const maxnodes=" + str(max_nodes) + "." + "\n\n" + "#const maxedges=" + str(max_edges) + "."

    final_text = nodes_to_avoid + "\n\n" + edges_to_avoid + "\n\n" + formulas + "\n\n" + dnfs + "\n\n" + clauses + "\n\n" + subspaces + "\n\n" + goals + "\n\n" + constants
    if filename != "":
        # Saving file
        with open(filename + ".asp", "w") as file:
            file.write(final_text.lower())
    return final_text.lower()


def is_included_in_subspace(subspace1: dict, subspace2: dict):
    """
    Test whether *subspace1* is contained in *subspace2*.

    **arguments**:
        * *subspace1*, *subspace2* (dicts): subspaces.
    **returns**:
        * Answer (bool): whether *subspace1* is contained in *subspace2*.
    **example**::
        >>> is_included_in_subspace({'v1': 0, 'v2': 1}, {'v2': 1})
        True
    """

    return all(x in subspace1 and subspace1[x] == subspace2[x] for x in subspace2.keys())


def select_trapspaces(tspaces, subspace: dict, use_attractors: bool = False, tsmin: List[dict] = None, complex_attractors: List[dict] = None):
    """
    Returns the trap spaces from *tspaces* that are contained in *subspace*.
    If *use_attractors* is True, it also returns the trap spaces from *tspaces* that contain only elements from *tsmin* and *complex_attractors* that are contained in *subspace*.
    It does not check that the elements of *tsmin* are minimal trap spaces or that the elements of *complex_attractors* are attractors.
    **arguments**:
        * *tspaces*: list of trap spaces.
        * *subspace*: subspace.
        * *use_attractors* (bool): indicates whether attractors are used in the selection of trap spaces or not. Default value: False.
        * *tsmin*: minimal trap spaces. Only used when *use_attractors* is True. Default value: [].
        * *complex_attractors*: list of complex attractors. A complex attractor is expected as a list of states (dicts). Only used when *use_attractors* is True. Default value: [].
    **returns**:
        * Selected (list): the trap spaces contained in *subspace* and, if *use_attractors* is True, also the trap spaces that contain only elements from *Tsmin* and *ComplexAttractors* that are contained in *subspace*.
    **example**::
        >>> selectTrapSpaces([{'v1': 0,'v2': 1}, {'v1': 0,'v3': 1}, {'v1': 0,'v2': 1,'v3': 1}], {'v2': 1})
        [{'v1': 0,'v2': 1}, {'v1': 0,'v2': 1,'v3': 1}]
    """

    if not tsmin:
        tsmin = []
    if not complex_attractors:
        complex_attractors = []

    # Trap spaces contained in *subspace*
    sel1 = [x for x in tspaces if is_included_in_subspace(x, subspace)]

    if not use_attractors:
        return sel1

    # Classify minimal trap spaces and complex attractors
    tsmin_accepted = [x for x in tsmin if is_included_in_subspace(x, subspace)]
    tsmin_discarded = [x for x in tsmin if x not in tsmin_accepted]
    cattr_accepted = [x for x in complex_attractors if all(is_included_in_subspace(y, subspace) for y in x)]
    cattr_discarded = [x for x in complex_attractors if x not in cattr_accepted]

    # If conditions cannot be matched
    if len(tsmin_accepted) + len(cattr_accepted) == 0:
        return sel1

    # If all trap spaces satisfy the condition
    if len(tsmin_discarded) + len(cattr_discarded) == 0:
        return tspaces

    tspaces_left = [x for x in tspaces if x not in sel1]
    sel2 = [ts for ts in tspaces_left
            if (any(is_included_in_subspace(x, ts) for x in tsmin_accepted) or
                any(is_included_in_subspace(y, ts) for x in cattr_accepted for y in x))
               and not any(is_included_in_subspace(y, ts) for y in tsmin_discarded)
               and not any(is_included_in_subspace(y, ts) for x in cattr_discarded for y in x)]

    return sel1 + sel2


def results_info(list_cs: List[dict]):
    """
    Returns a string stating the amount and size of the elements in *list_cs*.
    **returns**:
        * *text* (string): text stating the number and size of the elements in *CS*.
    **example**::
        >>> results_info([{'v1': 1}, {'v2':0, 'v3':1}])
    "2 control strategies, 1 of size 1, 1 of size 2"
    """

    text = str(len(list_cs)) + " control strategies"
    sizes = [len(x) for x in list_cs]
    cs_sizes = list({(el, sizes.count(el)) for el in sizes})
    cs_sizes.sort()
    for x in cs_sizes:
        text = text + ", " + str(x[1]) + " of size " + str(x[0])
        text = f"{len(list_cs)} control strategies, " + ", ".join(f"{x[1]} of size {x[0]}" for x in cs_sizes)
    return text



PRIME_JSON_FILE = "pyboolnet_primes.json"
cache = {}

def make_primes(bnfile, cachedir):
    cache_file = f"{cachedir}/{PRIME_JSON_FILE}"
    if cachedir and os.path.isfile(cache_file):
        try:
            with open(cache_file) as _f:
                primes = json.load(_f)
            main_logger.info("Loaded precomputed pyboolnet primes successfully")
            return primes
        except:
            main_logger.info("Loading precomputed pyboolnet primes fails")

    primes = time_check(bnet2primes)(bnfile)

    if cachedir:
        with open(cache_file, "w") as _f:
            json.dump(primes, _f)

    return primes


class PyBoolNet_ModelChecking:
    uses_cache = True
    bn_type = "bnet_file"

    @classmethod
    @time_check
    def run(self, bn: str, max_size:int, target:dict, exclude:list,
            expid:int, cachedir:str):
        if expid not in cache:
            cache[expid] = make_primes(bn, cachedir)

        return compute_control_strategies_with_model_checking(
                    primes=cache[expid],
                    avoid_nodes=exclude,
                    target=[target],
                    update=self.update)

    @staticmethod
    def free_experiment(expid):
        if expid in cache:
            del cache[expid]


@register_tool
class PyBoolNet_ModelChecking_SA(PyBoolNet_ModelChecking):
    name = "PBN[SA]"
    update = "synchronous"

@register_tool
class PyBoolNet_ModelChecking_ASA(PyBoolNet_ModelChecking):
    name = "PBN[ASA]"
    update = "asynchronous"


class PyBoolNet_Heuristic:
    uses_cache = True
    bn_type = "bnet_file"

    @classmethod
    @time_check
    def run(self, bn: str, max_size:int, target:dict, exclude:list,
            expid:int, cachedir:str):
        if expid not in cache:
            cache[expid] = make_primes(bn, cachedir)

        return run_control_problem(
                    primes=cache[expid],
                    avoid_nodes=exclude,
                    limit=max_size,
                    target=target,
                    control_type=self.control_type,
                    intervention_type="node")

    @staticmethod
    def free_experiment(expid):
        if expid in cache:
            del cache[expid]


@register_tool
class PyBoolNet_Percolation(PyBoolNet_Heuristic):
    name = "PBN[percolation]"
    control_type = "percolation"

@register_tool
class PyBoolNet_Trapspaces(PyBoolNet_Heuristic):
    name = "PBN[trap_spaces]"
    control_type = "trap_spaces"
