import os
import pandas as pd
from pathlib import Path
from collections import defaultdict
import re

# =========================
# CONFIG
# =========================
FOLDER_A = "/Users/victorferro/Documents/ufal/TCC_GIT_CPP/tcc_cpp/tcc_cpp/analises/ils_proccess_30_seg"
FOLDER_B = "/Users/victorferro/Documents/ufal/TCC_GIT_CPP/tcc_cpp/tcc_cpp/analises/ils_crude_30_seg"

OUTPUT_FILE = "per_instance_comparison.csv"

# =========================
# HELPERS
# =========================
def instance_name(filename):
    return filename.split("_seed_")[0]


def instance_sort_key(name):
    match = re.match(r"(\d+)(.*)", name)
    if match:
        return (int(match.group(1)), match.group(2))
    return (float("inf"), name)

def better_on_objective(avg_res_A, avg_res_B):
    if avg_res_B < avg_res_A:
        return "Crude"
    if avg_res_A < avg_res_B:
        return "Preprocessed"
    return "Tie"


def read_single_run(csv_path):
    df = pd.read_csv(csv_path, header=None)
    data = dict(zip(df[0], df[3]))

    objective = float(data["objective"])   # lower is better
    iterations = int(data["iterations"])   # higher is better

    return objective, iterations


def process_folder(folder):
    data = defaultdict(lambda: {
        "results": [],
        "iterations": []
    })

    for file in os.listdir(folder):
        if file.endswith(".csv"):
            inst = instance_name(file)
            obj, it = read_single_run(Path(folder) / file)
            data[inst]["results"].append(obj)
            data[inst]["iterations"].append(it)

    return data


# =========================
# PROCESS
# =========================
data_A = process_folder(FOLDER_A)
data_B = process_folder(FOLDER_B)

instances = sorted(set(data_A) & set(data_B), key=instance_sort_key)

rows = []

for inst in instances:
    # ---- Folder A ----
    avg_res_A = sum(data_A[inst]["results"]) / len(data_A[inst]["results"])
    best_res_A = min(data_A[inst]["results"])
    avg_it_A = sum(data_A[inst]["iterations"]) / len(data_A[inst]["iterations"])
    best_it_A = max(data_A[inst]["iterations"])

    # ---- Folder B ----
    avg_res_B = sum(data_B[inst]["results"]) / len(data_B[inst]["results"])
    best_res_B = min(data_B[inst]["results"])
    avg_it_B = sum(data_B[inst]["iterations"]) / len(data_B[inst]["iterations"])
    best_it_B = max(data_B[inst]["iterations"])

    # ---- Differences ----
    result_diff_abs = round(avg_res_B - avg_res_A,1)
    iteration_diff_abs = round(avg_it_A - avg_it_B,1)

    result_diff_pct = round((result_diff_abs / avg_res_A) * 100,2)
    iteration_diff_pct = round((iteration_diff_abs / avg_it_A) * 100,2)

    rows.append({
        # "instance": inst,
        # "avg_result_pre": avg_res_A,
        "best_result_pre": best_res_A,
        # "avg_iter_pre": avg_it_A,
        # "best_iter_pre": best_it_A,

        # "avg_result_crude": avg_res_B,
        # "best_result_crude": best_res_B,
        # "avg_iter_crude": avg_it_B,
        # "best_iter_crude": best_it_B,

        # "result_diff_abs": result_diff_abs,
        # "result_diff_%": result_diff_pct,

        # "iteration_diff_abs": iteration_diff_abs,
        # "iteration_diff_%": iteration_diff_pct,
        # "better_on_avg_objective": better_on_objective(avg_res_A, avg_res_B),
        # "better_on_best_objective": better_on_objective(best_res_A, best_res_B)
    })

# =========================
# OUTPUT
# =========================
df = pd.DataFrame(rows)
df.to_csv(OUTPUT_FILE, index=False)

print(f"Results written to: {OUTPUT_FILE}")
print(df)
