import os
import sys
import time
import csv
import warnings
import threading
import multiprocessing as mp
import pandas as pd
from tqdm.auto import tqdm
from micom import Community
from micom.solution import reset_solver

# Prevent thread thrashing across multiprocessing workers
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

LP_TIME_LIMIT = 2.0
STALL_WARN_SECONDS = 30

WORKER_COM = None
BM_S1_OBJ = None
BM_S2_OBJ = None
S2_GENE_MAP_GLOBAL = None
S2_ACTIVE_GENES_GLOBAL = None
HEARTBEAT = None
NAME1 = None
NAME2 = None


def classify_result(com):
    """Run com.optimize(); return 'optimal' or 'infeasible'."""
    sol = com.optimize()
    if sol is not None:
        return "optimal"

    try:
        reset_solver(com)
    except Exception:
        pass

    sol = com.optimize()
    if sol is not None:
        return "optimal"

    return "infeasible"


def init_worker(thresh_s1, bm_s1_base_id, name1, thresh_s2, bm_s2_base_id, name2, taxonomy_list, s2_gene_map=None, s2_active_genes=None, heartbeat=None):
    """Worker initialization function executed in each process pool child."""
    global WORKER_COM, BM_S1_OBJ, BM_S2_OBJ, S2_GENE_MAP_GLOBAL, S2_ACTIVE_GENES_GLOBAL, HEARTBEAT, NAME1, NAME2
    S2_GENE_MAP_GLOBAL = s2_gene_map
    S2_ACTIVE_GENES_GLOBAL = s2_active_genes
    HEARTBEAT = heartbeat
    NAME1 = name1
    NAME2 = name2
    warnings.filterwarnings("ignore")

    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    try:
        WORKER_COM = Community(pd.DataFrame(taxonomy_list))
        WORKER_COM.solver = 'gurobi'
        problem = WORKER_COM.solver.problem
        problem.Params.OutputFlag = 0
        problem.Params.LogToConsole = 0
        problem.Params.LogFile = ""
        problem.Params.Threads = 1
        problem.Params.TimeLimit = LP_TIME_LIMIT
        problem.Params.Method = 1
        problem.Params.Presolve = 0
        problem.Params.Crossover = 0

        BM_S1_OBJ = WORKER_COM.reactions.get_by_id(f"{bm_s1_base_id}__{NAME1}")
        BM_S2_OBJ = WORKER_COM.reactions.get_by_id(f"{bm_s2_base_id}__{NAME2}")

        BM_S1_OBJ.lower_bound = thresh_s1
        BM_S2_OBJ.lower_bound = thresh_s2
    except Exception as e:
        with open("worker_fatal_error.log", "a") as f:
            f.write(str(e) + "\n")
    finally:
        sys.stdout = original_stdout


def _tick():
    """Increment shared heartbeat counter."""
    if HEARTBEAT is not None:
        with HEARTBEAT.get_lock():
            HEARTBEAT.value += 1


def worker_single_lethal(args):
    """Worker function for single-gene lethal screens."""
    gene_id, rxn_ids = args
    global WORKER_COM

    orig_bounds = {}
    for r_id in rxn_ids:
        try:
            rxn = WORKER_COM.reactions.get_by_id(r_id)
            orig_bounds[r_id] = rxn.bounds
            rxn.bounds = (0.0, 0.0)
        except KeyError:
            pass

    try:
        result = classify_result(WORKER_COM)
        _tick()
    except Exception:
        result = "infeasible"
    finally:
        for r_id, b in orig_bounds.items():
            WORKER_COM.reactions.get_by_id(r_id).bounds = b

    if result == "infeasible":
        return gene_id
    return None


def worker_double_lethal(args):
    """Worker function for double-gene lethal screens with Fast-SL pruning."""
    gene_s1, rxn_ids_s1, thresh_s1, thresh_s2, is_s1_active = args
    global WORKER_COM, BM_S1_OBJ, BM_S2_OBJ, S2_GENE_MAP_GLOBAL, S2_ACTIVE_GENES_GLOBAL, NAME1, NAME2

    task_start = time.time()
    lethals_found = []

    orig_bounds_s1 = {}
    for r_id in rxn_ids_s1:
        try:
            rxn = WORKER_COM.reactions.get_by_id(r_id)
            orig_bounds_s1[r_id] = rxn.bounds
            rxn.bounds = (0.0, 0.0)
        except KeyError:
            pass

    try:
        for gene_s2, rxn_ids_s2 in S2_GENE_MAP_GLOBAL.items():
            if not is_s1_active and gene_s2 not in S2_ACTIVE_GENES_GLOBAL:
                continue

            orig_bounds_s2 = {}
            for r_id in rxn_ids_s2:
                try:
                    rxn = WORKER_COM.reactions.get_by_id(r_id)
                    orig_bounds_s2[r_id] = rxn.bounds
                    rxn.bounds = (0.0, 0.0)
                except KeyError:
                    pass

            try:
                result = classify_result(WORKER_COM)
                _tick()

                if result == "infeasible":
                    try:
                        BM_S2_OBJ.lower_bound = 0.0
                        s1_dead = classify_result(WORKER_COM) != "optimal"
                        _tick()

                        BM_S2_OBJ.lower_bound = thresh_s2
                        BM_S1_OBJ.lower_bound = 0.0
                        s2_dead = classify_result(WORKER_COM) != "optimal"
                        _tick()

                        if s1_dead and s2_dead:
                            lethals_found.append({"Target1": gene_s1, "Target2": gene_s2, "Lethality_Type": "Community Collapse"})
                        elif s1_dead:
                            lethals_found.append({"Target1": gene_s1, "Target2": gene_s2, "Lethality_Type": f"{NAME1} Death"})
                        else:
                            lethals_found.append({"Target1": gene_s1, "Target2": gene_s2, "Lethality_Type": f"{NAME2} Death"})
                    finally:
                        BM_S1_OBJ.lower_bound = thresh_s1
                        BM_S2_OBJ.lower_bound = thresh_s2
            except Exception:
                pass
            finally:
                for r_id, b in orig_bounds_s2.items():
                    WORKER_COM.reactions.get_by_id(r_id).bounds = b
    finally:
        for r_id, b in orig_bounds_s1.items():
            WORKER_COM.reactions.get_by_id(r_id).bounds = b

    elapsed = time.time() - task_start
    return (gene_s1, lethals_found, elapsed)


def start_watchdog(stop_event, heartbeat_val):
    """Silently tracks solves; only prints if there is a true stall."""
    last_count = 0
    last_check = time.time()
    stalled_since = None
    while not stop_event.is_set():
        time.sleep(5)
        now = time.time()
        with heartbeat_val.get_lock(): 
            current = heartbeat_val.value
        delta = current - last_count
        if delta == 0:
            if stalled_since is None: stalled_since = now
            stalled_secs = now - stalled_since
            if stalled_secs >= STALL_WARN_SECONDS:
                sys.stdout.write(f"\n   !! WATCHDOG: 0 LP solves in {stalled_secs:.0f}s. Stalled.\n")
                sys.stdout.flush()
        else:
            stalled_since = None
        last_count = current
        last_check = now


def run_pipeline(model1, name1, model2, name2, cores=None, death_fraction=0.01, output_csv=None, fraction=0.75, tolerance=1e-6):
    """Main entry point function to execute the Fast-SL workflow."""
    if cores is None:
        cores = max(1, mp.cpu_count() - 1)

    if output_csv is None:
        output_csv = f"{name1}_{name2}_Fastsl.csv"

    start_time = time.time()
    print(f"=== Fast-SL Community Screener ({name1} + {name2}) ===")
    print(f"Cores: {cores} | Death Threshold: {death_fraction}")

    taxonomy_list = [
        {"id": name1, "abundance": 0.5, "file": model1},
        {"id": name2, "abundance": 0.5, "file": model2}
    ]

    master_com = Community(pd.DataFrame(taxonomy_list))
    try:
        master_com.solver = 'gurobi'
        master_com.solver.problem.Params.OutputFlag = 0
        master_com.solver.problem.Params.Threads = 1
    except Exception:
        pass

    print("1. Calculating Wild-Type Cooperative Tradeoff (pFBA = True)...")
    wt_sol = master_com.cooperative_tradeoff(fraction=fraction, fluxes=True, pfba=True)
    wt_c = wt_sol.growth_rate
    wt_s1 = wt_sol.members.loc[name1, 'growth_rate']
    wt_s2 = wt_sol.members.loc[name2, 'growth_rate']

    print(f"   WT Growth (Community): {wt_c:.4f}")
    print(f"   WT Growth ({name1}): {wt_s1:.4f}")
    print(f"   WT Growth ({name2}): {wt_s2:.4f}\n")

    if wt_c < 1e-5 or wt_s1 < 1e-5 or wt_s2 < 1e-5:
        raise ValueError("Zero Wild-Type growth detected. Ensure viable media conditions.")

    thresh_s1 = death_fraction * wt_s1
    thresh_s2 = death_fraction * wt_s2

    s1_fluxes = wt_sol.fluxes.loc[name1].dropna()
    s2_fluxes = wt_sol.fluxes.loc[name2].dropna()

    bm_s1_id_base = s1_fluxes[(s1_fluxes - wt_s1).abs() < 1e-5].index[0].replace(f'__{name1}', '')
    bm_s2_id_base = s2_fluxes[(s2_fluxes - wt_s2).abs() < 1e-5].index[0].replace(f'__{name2}', '')

    active_s1_rxns = set(s1_fluxes[s1_fluxes.abs() > tolerance].index)
    active_s2_rxns = set(s2_fluxes[s2_fluxes.abs() > tolerance].index)

    print("2. Mapping Gene-Protein-Reaction (GPR) Associations...")
    s1_gene_map, s2_gene_map = {}, {}
    s1_active_genes, s2_active_genes = set(), set()

    for g in master_com.genes:
        if not g.reactions:
            continue

        saved = {r.id: r.bounds for r in g.reactions}
        g.knock_out()
        ko_rxns = [r.id for r in g.reactions if r.lower_bound == 0 and r.upper_bound == 0]
        g.functional = True
        for r_id, b in saved.items():
            master_com.reactions.get_by_id(r_id).bounds = b

        is_sp1 = any(f'__{name1}' in r.id for r in g.reactions)
        is_sp2 = any(f'__{name2}' in r.id for r in g.reactions)

        if is_sp1:
            if ko_rxns:
                s1_gene_map[g.id] = ko_rxns
            if any(r.replace(f'__{name1}', '') in active_s1_rxns for r in ko_rxns):
                s1_active_genes.add(g.id)
        elif is_sp2:
            if ko_rxns:
                s2_gene_map[g.id] = ko_rxns
            if any(r.replace(f'__{name2}', '') in active_s2_rxns for r in ko_rxns):
                s2_active_genes.add(g.id)

    print(f"   -> {name1} Genes: {len(s1_gene_map)} total ({len(s1_active_genes)} active)")
    print(f"   -> {name2} Genes: {len(s2_gene_map)} total ({len(s2_active_genes)} active)\n")

    ctx = mp.get_context('spawn')
    heartbeat_val = ctx.Value('L', 0)

    # --- SINGLE LETHALS ---
    s1_jsl, s2_jsl = [], []
    print(f"3. Running Single Lethal Screens on {cores} Cores...")

    with ctx.Pool(processes=cores, initializer=init_worker, initargs=(thresh_s1, bm_s1_id_base, name1, thresh_s2, bm_s2_id_base, name2, taxonomy_list, None, None, heartbeat_val), maxtasksperchild=50) as pool:
        args_s1 = [(g, rxns) for g, rxns in s1_gene_map.items() if g in s1_active_genes]
        for res in tqdm(pool.imap_unordered(worker_single_lethal, args_s1, chunksize=10), total=len(args_s1), desc=f"{name1} Singles"):
            if res:
                s1_jsl.append(res)

        args_s2 = [(g, rxns) for g, rxns in s2_gene_map.items() if g in s2_active_genes]
        for res in tqdm(pool.imap_unordered(worker_single_lethal, args_s2, chunksize=10), total=len(args_s2), desc=f"{name2} Singles"):
            if res:
                s2_jsl.append(res)

    print(f"   -> {name1} Single Lethals: {len(s1_jsl)}")
    print(f"   -> {name2} Single Lethals: {len(s2_jsl)}\n")

    # --- DOUBLE LETHALS ---
    s1_clean = {g: rxns for g, rxns in s1_gene_map.items() if g not in s1_jsl}
    s2_clean = {g: rxns for g, rxns in s2_gene_map.items() if g not in s2_jsl}

    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([f"{name1}_Gene", f"{name2}_Gene", "Lethality_Type"])

    total_pruned = (len(s1_clean) * len(s2_clean)) - (
        sum(1 for g in s1_clean if g not in s1_active_genes) *
        sum(1 for g in s2_clean if g not in s2_active_genes)
    )

    print(f"4. Running Double Lethal Screen (~{total_pruned} pairs evaluated)...")
    all_results = []

    with ctx.Pool(processes=cores, initializer=init_worker, initargs=(thresh_s1, bm_s1_id_base, name1, thresh_s2, bm_s2_id_base, name2, taxonomy_list, s2_clean, s2_active_genes, heartbeat_val), maxtasksperchild=50) as pool:
        args_double = [(g_s1, rxns, thresh_s1, thresh_s2, g_s1 in s1_active_genes) for g_s1, rxns in s1_clean.items()]
        
        stop_event = threading.Event()
        watchdog = threading.Thread(target=start_watchdog, args=(stop_event, heartbeat_val), daemon=True)
        watchdog.start()

        try:
            for g_s1_done, result_list, elapsed in tqdm(pool.imap_unordered(worker_double_lethal, args_double, chunksize=1), total=len(args_double), desc="Double Lethals"):
                if result_list:
                    all_results.extend(result_list)
                    with open(output_csv, mode='a', newline='') as file:
                        writer = csv.writer(file)
                        for res in result_list:
                            writer.writerow([res["Target1"], res["Target2"], res["Lethality_Type"]])
        finally:
            stop_event.set()
            watchdog.join(timeout=2)

    elapsed_total = (time.time() - start_time) / 60.0
    print(f"\n=== PIPELINE COMPLETE ({elapsed_total:.2f} mins) ===")
    print(f"Saved results to: {output_csv}\n")

    return pd.DataFrame(all_results)
