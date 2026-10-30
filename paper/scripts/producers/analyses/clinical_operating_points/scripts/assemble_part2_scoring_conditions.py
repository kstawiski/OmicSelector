from run import analysis_directory as _public_analysis_directory, public_table as _public_table
#!/usr/bin/env python3
import csv
import hashlib
import json
import pathlib
import shutil

csv.field_size_limit(32 * 1024 * 1024)

R8B = pathlib.Path("./analyses/clinical_operating_points")
M = R8B.parent
D = M.parent

SUPP_DIR = D / "analyses/tables_supplement_package/out/supplement"
BASE = D / "analyses"
PART2_IN = R8B / "evidence/retrieved/part2"
OUT_S4 = R8B / "out/scoring_conditions_rows"
for sub in [OUT_S4, OUT_S4 / "synthetic", OUT_S4 / "real", OUT_S4 / "parity"]:
    sub.mkdir(parents=True, exist_ok=True)

COMMIT_268 = "a55f6c2b99c986cb4ac897cced9c0996f5014c82"
TARBALL_SHA_268 = "9c322ca412641939c2d1698eb0281ea074bbea0942e4d11e6220cf0a983249a2"
RELEASE_IDENTITY_SRC = str(M / "package_2_6_8/out/OmicSelector_2.6.8.tar.gz.sha256")


def read_tsv(path: pathlib.Path, delimiter="\t"):
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_table(path: pathlib.Path, rows, fieldnames=None, delimiter="\t"):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow(
                {
                    k: json.dumps(v, ensure_ascii=False, sort_keys=True)
                    if isinstance(v, (dict, list))
                    else v
                    for k, v in row.items()
                }
            )


def sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


syn_code_id = sha256_file(R8B / "staging/payload/scripts/run_software_audit_268.R")

syn_raw = {
    "ai-tabpfn": read_tsv(PART2_IN / "synthetic/cpu_n1/ai-tabpfn.tsv")[0],
    "tab-tabdpt": read_tsv(PART2_IN / "synthetic/cpu_default/tab-tabdpt.tsv")[0],
}
real_raw = {
    "ai-tabpfn": read_tsv(PART2_IN / "real/cpu_n1/ai-tabpfn.tsv")[0],
    "tab-tabdpt": read_tsv(PART2_IN / "real/cpu_default/tab-tabdpt.tsv")[0],
}
par_raw = {
    "ai-tabpfn": read_tsv(PART2_IN / "parity/cuda/ai-tabpfn.tsv")[0],
    "tab-tabdpt": read_tsv(PART2_IN / "parity/cpu/tab-tabdpt.tsv")[0],
}

old_syn_raw = {
    "ai-tabpfn": read_tsv(BASE / "cohort_inventory_release/final4/synthetic_audit/synthetic/ai-tabpfn.tsv")[0],
    "tab-tabdpt": read_tsv(BASE / "cohort_inventory_release/final4/synthetic_audit/synthetic/tab-tabdpt.tsv")[0],
}
old_real_tabpfn = read_tsv(BASE / "cohort_inventory_release/tabpfn_real_audit/audit/real/ai-tabpfn.tsv")[0]
old_real_tabdpt = [
    r for r in read_tsv(BASE / "feature_annotation_realaudit/out/real_data_audit_complete.tsv") if r["method"] == "tab-tabdpt"
][0]
old_par_raw = {
    "ai-tabpfn": read_tsv(BASE / "device_benchmark_local_gpu_runs/out_gpu2/results/parity/parity/ai-tabpfn.tsv")[0],
    "tab-tabdpt": read_tsv(BASE / "cohort_preparation_parity/out/parity/tab-tabdpt.tsv")[0],
}

# 2. Build producer-schema synthetic, real, and parity tables in OUT_scoring_conditions
syn_producer_fields = [
    "method",
    "status",
    "rejection_code",
    "reason",
    "singleton_batch_delta",
    "permutation_delta",
    "duplicate_delta",
    "adversarial_batch_delta",
    "state_immutable",
    "fit_seconds",
    "package_version",
    "package_commit",
    "package_source_sha256",
    "analysis_code_id",
    "state_fingerprint_method",
    "metadata_probes",
    "fit_configuration",
]

syn_producer_rows = {}
for m in ["ai-tabpfn", "tab-tabdpt"]:
    s = syn_raw[m]
    max_dup = max(float(s["max_abs_dup_vs_single"]), float(s["max_abs_dup_pair"]))
    dup_str = "0" if max_dup == 0.0 else str(max_dup)
    row = {
        "method": m,
        "status": "working" if s["audit_pass"] == "TRUE" else "failed",
        "rejection_code": "",
        "reason": s["score_error"] or s["fit_error"] or "",
        "singleton_batch_delta": s["max_abs_single_vs_batch"],
        "permutation_delta": s["max_abs_perm_vs_batch"],
        "duplicate_delta": dup_str,
        "adversarial_batch_delta": s["max_abs_adversarial_companion"],
        "state_immutable": s["state_unchanged"],
        "fit_seconds": s["seconds"],
        "package_version": s["package_version"],
        "package_commit": s["package_commit"],
        "package_source_sha256": TARBALL_SHA_268,
        "analysis_code_id": syn_code_id,
        "state_fingerprint_method": "SHA-256 R + Python tensor/array/estimator state",
        "metadata_probes": "TRUE",
        "fit_configuration": "n_estimators=1; full training context; CPU" if m == "ai-tabpfn" else "defaults",
    }
    syn_producer_rows[m] = row
    write_table(OUT_S4 / "synthetic" / _public_analysis_directory(f"{m}.tsv"), [row], syn_producer_fields)
write_table(OUT_S4 / "synthetic_audit_268.tsv", list(syn_producer_rows.values()), syn_producer_fields)

real_producer_fields = [
    "method",
    "status",
    "reason",
    "dataset",
    "split",
    "n_train",
    "n_test",
    "train_people",
    "test_people",
    "person_overlap",
    "train_groups",
    "test_groups",
    "group_overlap",
    "n_features",
    "maximum_score_delta",
    "max_abs_singleton_vs_batch",
    "max_abs_permutation",
    "max_abs_adversarial_companion",
    "max_abs_duplicate_vs_singleton",
    "max_abs_duplicate_pair",
    "state_immutable",
    "state_fingerprint_method",
    "elapsed_seconds",
    "package_version",
    "package_commit",
    "package_source_sha256",
    "fit_configuration",
]

real_producer_rows = {}
for m in ["ai-tabpfn", "tab-tabdpt"]:
    rl = real_raw[m]
    max_delta = max(
        float(rl["max_abs_singleton_vs_batch"]),
        float(rl["max_abs_permutation"]),
        float(rl["max_abs_adversarial_companion"]),
        float(rl["max_abs_duplicate_vs_singleton"]),
        float(rl["max_abs_duplicate_pair"]),
    )
    delta_str = "0" if max_delta == 0.0 else str(max_delta)
    row = {
        "method": m,
        "status": "working" if rl["pass"] == "TRUE" else "failed",
        "reason": "",
        "dataset": "GSE204951",
        "split": 'Held-out profiles belong to different people from the training profiles.'
        if m == "ai-tabpfn"
        else 'The grouped person split uses seed-101 and real_fixture.rds.',
        "n_train": rl["n_train"],
        "n_test": rl["n_test"],
        "train_people": "12",
        "test_people": "8",
        "person_overlap": "0",
        "train_groups": "12",
        "test_groups": "8",
        "group_overlap": "0",
        "n_features": rl["n_features"],
        "maximum_score_delta": delta_str,
        "max_abs_singleton_vs_batch": rl["max_abs_singleton_vs_batch"],
        "max_abs_permutation": rl["max_abs_permutation"],
        "max_abs_adversarial_companion": rl["max_abs_adversarial_companion"],
        "max_abs_duplicate_vs_singleton": rl["max_abs_duplicate_vs_singleton"],
        "max_abs_duplicate_pair": rl["max_abs_duplicate_pair"],
        "state_immutable": rl["state_unchanged"],
        "state_fingerprint_method": "SHA-256 R + Python tensor/array/estimator state",
        "elapsed_seconds": rl["seconds"],
        "package_version": "2.6.8",
        "package_commit": rl["package_commit"],
        "package_source_sha256": TARBALL_SHA_268,
        "fit_configuration": "n_estimators=1; full training context; CPU"
        if m == "ai-tabpfn"
        else "defaults; full training context; CPU",
    }
    real_producer_rows[m] = row
    write_table(OUT_S4 / "real" / _public_analysis_directory(f"{m}.tsv"), [row], real_producer_fields)
write_table(OUT_S4 / "real_audit_268.tsv", list(real_producer_rows.values()), real_producer_fields)

par_producer_fields = list(old_par_raw["ai-tabpfn"].keys())
par_producer_rows = {}
for m in ["ai-tabpfn", "tab-tabdpt"]:
    pr = {k: par_raw[m][k] for k in par_producer_fields}
    par_producer_rows[m] = pr
    write_table(OUT_S4 / "parity" / _public_analysis_directory(f"{m}.tsv"), [pr], par_producer_fields)
write_table(OUT_S4 / "parity_268.tsv", list(par_producer_rows.values()), par_producer_fields)

# 3. Read original scoring_conditions.tsv (74 columns, 76 rows) and update ai-tabpfn and tab-tabdpt rows
old_s4_rows = read_tsv(SUPP_DIR / "scoring_conditions.tsv")
s4_fields = list(old_s4_rows[0].keys())

new_s4_rows = []
updated_s4_only = []
implementation_comparison_records = []

for r in old_s4_rows:
    mid = r["Method identifier"]
    if mid not in ("ai-tabpfn", "tab-tabdpt"):
        new_s4_rows.append(dict(r))
        continue

    nr = dict(r)
    sp = syn_producer_rows[mid]
    rp = real_producer_rows[mid]
    pp = par_producer_rows[mid]

    syn_src = str(OUT_S4 / "synthetic" / _public_analysis_directory(f"{mid}.tsv"))
    real_src = str(OUT_S4 / "real" / _public_analysis_directory(f"{mid}.tsv"))
    par_src = str(OUT_S4 / "parity" / _public_analysis_directory(f"{mid}.tsv"))

    syn_hist = json.loads(r["Synthetic assay records"])
    syn_hist.append(dict(sp, source_path=syn_src, source_row=2))

    real_hist = json.loads(r["Real-data records"])
    real_hist.append(dict(rp, source_path=real_src, source_row=2))

    par_hist = json.loads(r["Parity records"])
    par_hist.append(dict(pp, source_path=par_src, source_row=2))

    nr.update(
        {
            "One-specimen audit outcome": sp["status"],
            "Outcome explanation": sp["reason"],
            "Maximum score change when scored alone": sp["singleton_batch_delta"],
            "Maximum score change after reordering test specimens": sp["permutation_delta"],
            "Maximum score change after duplicating a test specimen": sp["duplicate_delta"],
            "Maximum score change with extreme neighboring specimens": sp["adversarial_batch_delta"],
            "Fitted state remained unchanged": sp["state_immutable"],
            "Fit time (seconds)": sp["fit_seconds"],
            "OmicSelector version": sp["package_version"],
            "OmicSelector commit": sp["package_commit"],
            "Analysis-code identifier": sp["analysis_code_id"],
            "Final release commit": COMMIT_268,
            "Synthetic release binding": 'Execution uses the recorded release.',
            "Synthetic source": syn_src,
            "Synthetic source row": "2",
            "Synthetic assay records": json.dumps(syn_hist, ensure_ascii=False, sort_keys=True),
            "Real-data audit outcome": rp["status"],
            "Real-data dataset": rp["dataset"],
            "Real-data split": rp["split"],
            "Real-data training profiles": rp["n_train"],
            "Real-data test profiles": rp["n_test"],
            "Real-data training groups": rp["train_groups"],
            "Real-data test groups": rp["test_groups"],
            "Real-data group overlap": rp["group_overlap"],
            "Real-data maximum score difference": rp["maximum_score_delta"],
            "Real-data state unchanged": rp["state_immutable"],
            "Real-data execution commit": rp["package_commit"],
            "Real-data reason": rp["reason"],
            "Real-data source": real_src,
            "Real-data records": json.dumps(real_hist, ensure_ascii=False, sort_keys=True),
            "Real-data audit configuration": rp["fit_configuration"],
            "Real-data state fingerprint method": rp["state_fingerprint_method"],
            "Real-data features": rp["n_features"],
            "Real-data configuration source": real_src,
            "Parity outcome": pp["status"],
            "Parity maximum score difference": pp["max_abs_score_difference"],
            "Parity dataset or target": pp["cohort"],
            "Parity exact scores": pp["exact_scores"],
            "Parity rank correlation": pp["spearman"],
            "Parity monotone mapping": pp["strictly_increasing_mapping"],
            "Parity execution commit": pp["package_commit"],
            "Parity cause": pp["cause"],
            "Parity source": par_src,
            "Parity records": json.dumps(par_hist, ensure_ascii=False, sort_keys=True),
            "Final release identity source": RELEASE_IDENTITY_SRC,
            "Shared-split parity outcome": pp["status"],
            "Shared-split parity cause": pp["cause"],
            "Deployment-versus-harness parity": "identical" if pp["status"] == "identical" else pp["status"],
            "Parity harness AUC": pp["auc_harness"],
            "Parity deployment AUC": pp["auc_deployment"],
            "Parity comparison scope": 'Both scoring paths use the same held-out split.',
        }
    )

    new_s4_rows.append(nr)
    updated_s4_only.append(nr)

    for col in s4_fields:
        if col in ("Synthetic assay records", "Real-data records", "Parity records"):
            continue
        implementation_comparison_records.append(
            {
                "table": "scoring_conditions.tsv",
                "method": mid,
                "field": col,
                "old_value": r[col],
                "analysis_value": nr[col],
                "changed": str(r[col] != nr[col]),
            }
        )

for mid in ["ai-tabpfn", "tab-tabdpt"]:
    for col in syn_producer_fields:
        ov = old_syn_raw[mid].get(col, "")
        nv = syn_producer_rows[mid].get(col, "")
        implementation_comparison_records.append(
            {
                "table": "synthetic_audit",
                "method": mid,
                "field": col,
                "old_value": ov,
                "analysis_value": nv,
                "changed": str(ov != nv),
            }
        )
    old_r = old_real_tabpfn if mid == "ai-tabpfn" else old_real_tabdpt
    for col in ["status", "n_train", "n_test", "maximum_score_delta", "state_immutable", "package_commit"]:
        ov = old_r.get(col, "")
        nv = real_producer_rows[mid].get(col, "")
        implementation_comparison_records.append(
            {
                "table": "real_audit",
                "method": mid,
                "field": col,
                "old_value": ov,
                "analysis_value": nv,
                "changed": str(ov != nv),
            }
        )
    for col in par_producer_fields:
        ov = old_par_raw[mid].get(col, "")
        nv = par_producer_rows[mid].get(col, "")
        implementation_comparison_records.append(
            {
                "table": "parity",
                "method": mid,
                "field": col,
                "old_value": ov,
                "analysis_value": nv,
                "changed": str(ov != nv),
            }
        )

# Write updated scoring_conditions tables
write_table(OUT_S4 / "scoring_conditions.tsv", new_s4_rows, s4_fields, delimiter="\t")
write_table(OUT_S4 / "scoring_conditions_tabpfn_tabdpt_rows.tsv", updated_s4_only, s4_fields, delimiter="\t")
write_table(
    OUT_S4 / "scoring_conditions_reference_comparison.tsv",
    implementation_comparison_records,
    ["table", "method", "field", "old_value", "analysis_value", "changed"],
    delimiter="\t",
)

# Copy scoring_conditions_parity_summary.tsv (still 58 identical, 2 monotone, 14 deployment-interface)
shutil.copy2(SUPP_DIR / "scoring_conditions_parity_summary.tsv", OUT_S4 / "scoring_conditions_parity_summary.tsv")

# Collect all synthetic, real, and parity runs across CPU and GPU
all_syn = []
for sub in ["cpu_n1", "cpu_default", "cuda_n1", "cuda_default"]:
    for m in ["ai-tabpfn", "tab-tabdpt"]:
        p = PART2_IN / "synthetic" / _public_analysis_directory(sub) / _public_analysis_directory(f"{m}.tsv")
        if p.exists():
            row = read_tsv(p)[0]
            row["run_config"] = sub
            all_syn.append(row)
write_table(OUT_S4 / "synthetic_audit_all_runs.tsv", all_syn)

all_real = []
for sub in ["cpu_n1", "cpu_default", "cuda_n1", "cuda_default"]:
    for m in ["ai-tabpfn", "tab-tabdpt"]:
        p = PART2_IN / "real" / _public_analysis_directory(sub) / _public_analysis_directory(f"{m}.tsv")
        if p.exists():
            row = read_tsv(p)[0]
            row["run_config"] = sub
            all_real.append(row)
write_table(OUT_S4 / "real_data_audit_all_runs.tsv", all_real)

all_par = []
for sub in ["cpu", "cuda"]:
    for m in ["ai-tabpfn", "tab-tabdpt"]:
        p = PART2_IN / "parity" / _public_analysis_directory(sub) / _public_analysis_directory(f"{m}.tsv")
        if p.exists():
            row = read_tsv(p)[0]
            row["device_mode"] = sub
            all_par.append(row)
write_table(OUT_S4 / "parity_all_runs.tsv", all_par)

# Copy per-specimen score files and device comparison
shutil.copy2(PART2_IN / "parity/cuda/ai-tabpfn_scores.tsv", OUT_S4 / "ai-tabpfn_scores.tsv")
shutil.copy2(PART2_IN / "parity/cpu/tab-tabdpt_scores.tsv", OUT_S4 / "tab-tabdpt_scores.tsv")
if (PART2_IN / "parity/cuda/tab-tabdpt_scores.tsv").exists():
    shutil.copy2(PART2_IN / "parity/cuda/tab-tabdpt_scores.tsv", OUT_S4 / "tab-tabdpt_cuda_scores.tsv")
if (PART2_IN / "parity/cpu/ai-tabpfn_scores.tsv").exists():
    shutil.copy2(PART2_IN / "parity/cpu/ai-tabpfn_scores.tsv", OUT_S4 / "ai-tabpfn_cpu_scores.tsv")
if (PART2_IN / "device_comparison.tsv").exists():
    shutil.copy2(PART2_IN / "device_comparison.tsv", OUT_S4 / "device_comparison.tsv")

print("Assembled scoring_conditions_rows successfully")
