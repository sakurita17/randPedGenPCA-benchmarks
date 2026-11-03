"""Build auditable comparison tables from downloaded Eddie CSVs (stdlib only)."""
import argparse
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
from statistics import median


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def metric(row, field, divisor=1):
    value = number(row.get(field))
    return None if value is None else value / divisor


def key(row):
    return tuple(row.get(k, "") for k in ("n", "centredH", "q", "l", "seed"))


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or ["nota"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def markdown(rows, columns):
    def fmt(value, field):
        if value is None or value == "":
            return "—"
        if isinstance(value, float):
            if field.startswith("correlacion"):
                return f"{value:.9f}"
            return f"{value:.3g}"
        return str(value).replace("|", "/")
    return "\n".join([
        "| " + " | ".join(title for _, title in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
        *["| " + " | ".join(fmt(row.get(k), k) for k, _ in columns) + " |" for row in rows]
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("eddie_results"))
    parser.add_argument("--output", type=Path, default=Path("reports/eddie_comparison"))
    parser.add_argument("--language", choices=("en", "es", "both"), default="en")
    args = parser.parse_args()
    if not args.input.is_dir():
        parser.error("input directory does not exist")
    args.output.mkdir(parents=True, exist_ok=True)
    data_dir = args.output / "data"
    data_dir.mkdir(exist_ok=True)
    scale, geno, detailed, references, failures, inventory, all_rows = ([] for _ in range(7))
    provenance = []
    for path in sorted(args.input.glob("*/results.csv")):
        run = path.parent.name
        config_path = path.parent / "config.json"
        config = json.loads(config_path.read_text())
        rows = list(csv.DictReader(path.open()))
        provenance.append(dict(run=run, file=str(path),
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                               config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest()))
        experiment = config.get("experiment", "large")
        if experiment == "large" and max(config["sizes"]) < 10000:
            experiment = "smoke"
        hashes = {Path(k).name: v for k, v in config.get("source_sha256", {}).items()
                  if Path(k).name in ("utils.R", "rhpca.R", "hutchH.R")}
        source_id = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()[:12]
        for row in rows:
            all_rows.append(dict(run=run, experiment=experiment, **row))
            if row.get("status") != "ok":
                failures.append(dict(run=run, experiment=experiment, **row))
        pcas = [row for row in rows if row["stage"] == "pca" and row["status"] == "ok"]
        diagnostics = {key(row): row for row in rows
                       if row["stage"] in ("diagnose", "compare_exact") and row["status"] == "ok"}
        # Genotyping levels share n/center/q/l/seed; keep ng in their join key.
        if experiment == "genotyped":
            diagnostics = {(key(row), row["n_genotyped"]): row for row in rows
                           if row["stage"] == "diagnose" and row["status"] == "ok"}
        if experiment == "exact":
            expected = len(set(config["sizes"])) * 2 * 2 * config["replicates"]
        elif experiment == "genotyped":
            expected = len(set(config["genotyped"])) * 2 * config["replicates"]
        else:
            centers = 2 if config["centering"] == "both" else 1
            expected = len(set(config["sizes"])) * centers * (2 if config["compare"] else 1) * config["replicates"]
        inventory.append(dict(run=run, experimento=experiment, pca_esperadas=expected,
                              pca_ok=len(pcas), diagnosticos_ok=len(diagnostics),
                              etapas_con_fallo=sum(r["status"] != "ok" for r in rows),
                              cobertura="completa" if len(pcas) == expected and len(diagnostics) == expected
                              else "incompleta", fuentes_id=source_id))
        for pca in pcas:
            diag = diagnostics.get((key(pca), pca["n_genotyped"]) if experiment == "genotyped" else key(pca), {})
            base = dict(run=run, fuentes_id=source_id, n=int(pca["n"]),
                        genotipados=int(pca["n_genotyped"]), marcadores=int(pca["markers"]),
                        operador="CHC" if pca["centredH"] == "TRUE" else "H",
                        q=int(pca["q"]), l=int(pca["l"]), semilla=int(pca["seed"]),
                        tiempo_pca_s=metric(pca, "pca_elapsed_sec"),
                        tiempo_pca_min=metric(pca, "pca_elapsed_sec", 60),
                        tiempo_pca_h=metric(pca, "pca_elapsed_sec", 3600),
                        proceso_pca_s=metric(pca, "process_wall_sec"),
                        rss_proceso_mib=metric(pca, "peak_process_rss_mib"),
                        rss_proceso_gib=metric(pca, "peak_process_rss_mib", 1024),
                        diagnostico_s=metric(diag, "diagnostic_elapsed_sec"),
                        residuo_max_relativo=metric(diag, "max_relative_eigenpair_residual"),
                        cambio_eigenvalores_vs_base=metric(diag, "eigenvalue_relative_change_vs_base"),
                        cambio_subespacio_vs_base=metric(diag, "normalised_subspace_change_vs_base"),
                        error_ortogonalidad=metric(pca, "orth_fro_error"),
                        score_check_status="aligned_U_d" if pca.get("score_scaling") == "U*d"
                          else "historical_check_requires_review",
                        advertencias=pca.get("warnings", ""),
                        diagnostico_estado=diag.get("status", "ausente"),
                        archivo_resultados=str(path))
            if experiment == "large":
                scale.append(base)
            elif experiment == "genotyped":
                base["porcentaje_genotipado"] = 100 * base["genotipados"] / base["n"]
                geno.append(base)
            elif experiment == "exact":
                # Old reports compared U*d scores with an exact U*sqrt(d)
                # reference. Keep their raw metric for audit, not as valid error.
                score_reference_valid = diag.get("reference_score_scaling") == "U*d"
                raw_score_error = metric(diag, "aligned_scores_relative_error")
                base.update(error_eigenvalores_pct=metric(diag, "eigenvalue_max_rel_error", .01),
                            error_subespacio=metric(diag, "normalised_subspace_error"),
                            error_scores_alineados=raw_score_error if score_reference_valid else None,
                            raw_historical_aligned_scores_error=raw_score_error,
                            score_comparison_status="aligned_U_d" if score_reference_valid
                              else "historical_reference_scale_not_aligned",
                            correlacion_pc1=metric(diag, "pc1_abs_cor"),
                            correlacion_pc2=metric(diag, "pc2_abs_cor"),
                            error_traza_pct=metric(diag, "trace_relative_error", .01),
                            error_proporcion_variacion=metric(diag, "variance_proportion_max_abs_error"))
                detailed.append(base)
        for ref in [r for r in rows if r["stage"] == "exact" and r["status"] == "ok"]:
            references.append(dict(run=run, n=int(ref["n"]),
                                   genotipados=int(ref["n_genotyped"]),
                                   operador="CHC" if ref["centredH"] == "TRUE" else "H",
                                   construir_H_s=metric(ref, "getH_sec"),
                                   eigen_exacto_s=metric(ref, "exact_eigen_sec"),
                                   H_mas_eigen_s=metric(ref, "getH_sec") + metric(ref, "exact_eigen_sec"),
                                   proceso_referencia_s=metric(ref, "process_wall_sec"),
                                   rss_referencia_mib=metric(ref, "peak_process_rss_mib"),
                                   eigenvalor_minimo=metric(ref, "exact_min_eigenvalue"),
                                   separacion_relativa_PC5_PC6=metric(ref, "relative_eigengap_at_rho")))
    if not inventory:
        parser.error("no run results found")
    scale.sort(key=lambda r: (r["n"], r["operador"] == "CHC", r["q"], r["semilla"], r["run"]))
    geno.sort(key=lambda r: (r["genotipados"], r["operador"] == "CHC", r["semilla"]))
    groups = defaultdict(list)
    for row in detailed:
        groups[(row["run"], row["n"], row["operador"], row["q"], row["l"])].append(row)
    exact = []
    for (run, n, operator, q, sketch), group in sorted(groups.items()):
        reference = next(r for r in references if r["run"] == run and r["n"] == n and r["operador"] == operator)
        summary = dict(run=run, n=n, genotipados=group[0]["genotipados"], operador=operator,
                       q=q, l=sketch, semillas=len({r["semilla"] for r in group}),
                       tiempo_pca_mediana_s=median(r["tiempo_pca_s"] for r in group),
                       tiempo_pca_min_s=min(r["tiempo_pca_s"] for r in group),
                       tiempo_pca_max_s=max(r["tiempo_pca_s"] for r in group),
                       referencia_H_mas_eigen_s=reference["H_mas_eigen_s"],
                       rss_pca_max_mib=max(r["rss_proceso_mib"] for r in group),
                       error_eigenvalores_peor_pct=max(r["error_eigenvalores_pct"] for r in group),
                       error_subespacio_peor=max(r["error_subespacio"] for r in group),
                       error_scores_alineados_peor=max(r["error_scores_alineados"] for r in group)
                         if all(r["error_scores_alineados"] is not None for r in group) else None,
                       score_comparison_status="aligned_U_d" if all(
                         r["error_scores_alineados"] is not None for r in group)
                         else "historical_reference_scale_not_aligned",
                       residuo_peor=max(r["residuo_max_relativo"] for r in group),
                       correlacion_pc1_min=min(r["correlacion_pc1"] for r in group),
                       correlacion_pc2_min=min(r["correlacion_pc2"] for r in group),
                       error_traza_abs_peor_pct=max(abs(r["error_traza_pct"]) for r in group),
                       error_proporcion_variacion_peor=max(r["error_proporcion_variacion"] for r in group))
        exact.append(summary)
    # Preserve successful repetitions rather than silently selecting a preferred run.
    successful_keys = defaultdict(list)
    for row in scale:
        successful_keys[tuple(row[k] for k in ("fuentes_id", "n", "genotipados", "marcadores", "operador", "q", "l", "semilla"))].append(row)
    duplicates = [row for group in successful_keys.values() if len(group) > 1 for row in group]
    # Audit nested animal IDs using the transferred membership CSVs.
    nested_checks = []
    for directory in sorted(args.input.glob("followup_genotyped*")):
        previous = set()
        files = sorted(directory.glob("data_g*.rds_ids.csv"), key=lambda p: int(p.name.split("g", 1)[1].split(".")[0]))
        for path in files:
            ids = [r["id"] for r in csv.DictReader(path.open())]
            expected = int(path.name.split("g", 1)[1].split(".")[0])
            passed = len(ids) == expected and len(set(ids)) == expected and previous <= set(ids)
            nested_checks.append(dict(run=directory.name, genotipados=expected,
                                      numero_ids=len(ids), anidacion_correcta=passed))
            if not passed:
                raise ValueError(f"Invalid nested IDs: {path}")
            previous = set(ids)
    tables = {"escalabilidad": scale, "precision_resumen": exact,
              "genotipados": geno, "precision_semillas": detailed,
              "referencia_exacta": references, "fallos": failures,
              "inventario": inventory, "duplicados": duplicates,
              "auditoria_anidacion": nested_checks, "etapas_originales": all_rows}
    for name, rows in tables.items():
        write_csv(data_dir / (name + ".csv"), rows)
    (data_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))
    report = ["# Comparación de resultados de Eddie", "",
              "Tablas generadas desde los CSV descargados. Las cifras describen estos experimentos; no prueban una ley general de escalabilidad.", "",
              "## Escalabilidad", "",
              markdown(scale, [("n", "Individuos"), ("operador", "Operador"), ("q", "q"), ("l", "l"),
                               ("tiempo_pca_h", "PCA (h)"), ("rss_proceso_gib", "RSS (GiB)"),
                               ("residuo_max_relativo", "Residuo relativo")]), "",
              "Solo PCA completadas; las repeticiones exitosas se conservan y se auditan en duplicados.csv. Las pruebas pequeñas de circuito se excluyen de este cuadro.", "",
              "## Precisión frente a referencia exacta", "",
              markdown(exact, [("n", "Individuos"), ("operador", "Operador"), ("q", "q"), ("l", "l"),
                               ("tiempo_pca_mediana_s", "PCA mediana (s)"),
                               ("error_eigenvalores_peor_pct", "Peor error eigenvalores (%)"),
                               ("error_subespacio_peor", "Peor error subespacio"),
                               ("correlacion_pc2_min", "Mín. correlación PC2")]), "",
              "Por configuración: mediana de tiempo y peor error de tres semillas de PCA sobre el mismo dataset. No son tres simulaciones independientes. El error de eigenvalores es el máximo relativo entre los cinco componentes, expresado como porcentaje. El error de subespacio es adimensional y no se expresa como porcentaje.", "",
              "La referencia diagonaliza exactamente H o CHC, no la matriz de genotipos. H_mas_eigen_s incluye construir H y la diagonalización completa; rhpca calcula cinco componentes e incluye Hutch++. No debe interpretarse su cociente como una comparación universal de algoritmos equivalentes.", "",
              "## Pedigrí fijo y más genotipados", "",
              markdown(geno, [("n", "Pedigrí"), ("genotipados", "Genotipados"), ("operador", "Operador"),
                              ("tiempo_pca_min", "PCA (min)"), ("rss_proceso_mib", "RSS (MiB)"),
                              ("residuo_max_relativo", "Residuo relativo")]), "",
              "q=4, l=20 y una semilla. La disminución observada de tiempo no demuestra que aumentar genotipados acelere generalmente el método: hacen falta repeticiones y control del entorno. El cambio de genotipados también cambia el operador; el residuo no es error frente a una PCA exacta.", "",
              "Los archivos de IDs verifican conjuntos anidados. La igualdad de pedigrí/genotipos fue comprobada por el generador; no se ha reverificado aquí porque los datasets .rds no se descargaron.", "",
              "## Cobertura de ejecuciones", "",
              markdown(inventory, [("run", "Ejecución"), ("pca_esperadas", "PCA previstas"),
                                   ("pca_ok", "PCA ok"), ("diagnosticos_ok", "Diagnósticos ok"),
                                   ("cobertura", "Cobertura")]), "",
              "No hay resultados descargados de CHC con 100.000 individuos. No se puede determinar solo con estos CSV si esas dos configuraciones se enviaron o qué ocurrió con ellas.", "",
              "El timeout explícito se conserva en fallos.csv y no se usa como tiempo de PCA completada. Una terminación del gestor puede no dejar fila de fallo en results.csv; se necesita qacct/log para conocer la causa. En la conversación se confirmó el corte de 48 horas de large_20261029_122040_658611.", "",
              "## Unidades y límites", "",
              "- tiempo_pca_s/min/h: tiempo interno de cálculo, incluida factorización del pedigrí y Hutch++, excluidos carga, guardado y diagnóstico.",
              "- proceso_pca_s: duración del proceso completo.",
              "- rss_proceso_mib/gib: máximo RSS del proceso completo; no memoria incremental exclusiva de PCA.",
              "- residuo: max_j ||K u_j - d_j u_j|| / ||K u_j||. No es el error contra eigenvalores exactos.",
              "- Los porcentajes de variación dependen de la traza estimada por Hutch++; eigenvalores precisos no garantizan porcentajes igual de precisos.",
              "- Todos los cuadros conservan ejecución y fuente para permitir auditoría. No hay información suficiente para atribuir diferencias de tiempo a una causa específica del nodo.", "",
              "## Regenerar en el Mac", "", "Desde la raíz del repositorio:", "", "```sh",
              "python3 validation/summarise_eddie.py",
              "Rscript validation/export_eddie_tables.R", "```", "",
              "El primer comando utiliza solamente Python estándar; la exportación Excel y gráficos requiere los paquetes R openxlsx y ggplot2, disponibles en este equipo.", ""]
    if args.language in ("es", "both"):
        (data_dir / "comparacion.md").write_text("\n".join(report))
    english = ["# Comparison of Eddie benchmark results", "",
               "Tables generated from the downloaded CSV files. These measurements describe these experiments; they do not establish a general scaling law.", "",
               "## Scalability", "",
               markdown(scale, [("n", "Individuals"), ("operador", "Operator"), ("q", "q"), ("l", "l"),
                                ("tiempo_pca_h", "PCA (h)"), ("rss_proceso_gib", "RSS (GiB)"),
                                ("residuo_max_relativo", "Relative residual")]), "",
               "Only successfully completed PCA runs are included. Successful repeated measurements are retained and audited separately. Small pipeline smoke tests are excluded from this table.", "",
               "## Accuracy against an exact reference", "",
               markdown(exact, [("n", "Individuals"), ("operador", "Operator"), ("q", "q"), ("l", "l"),
                                ("tiempo_pca_mediana_s", "Median PCA time (s)"),
                                ("error_eigenvalores_peor_pct", "Worst eigenvalue error (%)"),
                                ("error_subespacio_peor", "Worst subspace error"),
                                ("correlacion_pc2_min", "Minimum PC2 correlation")]), "",
               "Each configuration reports the median runtime and worst error across three PCA seeds applied to the same dataset. These are not three independent population simulations. Eigenvalue error is the maximum relative error among the five components, expressed as a percentage. Subspace error is dimensionless and is not expressed as a percentage. Correlations are absolute values, allowing arbitrary eigenvector signs; individual-PC correlations can be misleading when eigenvalues are nearly repeated.", "",
               "Score convention correction: the algorithm returns X = U diag(d). The original validation compared these scores with an exact U diag(sqrt(d)) reference. Historical aligned-score errors are therefore marked as not aligned and excluded from the valid error summary; raw values remain available for audit. Historical score_scaling_max_error checks also used sqrt(d). The validator has been corrected for future runs. Eigenvalue errors, subspace errors, residuals, absolute PC correlations and measured PCA times are unaffected by this score-scale correction. Recomputing coordinate-error diagnostics requires the saved RDS files, not a rerun of the expensive PCA.", "",
               "The exact reference diagonalizes the same H or CHC operator, not the genotype matrix. Reference runtime in the workbook includes building H and its full eigendecomposition; rhpca computes five components and includes Hutch++ trace estimation. A runtime ratio should not be interpreted as a universal comparison of equivalent algorithmic workloads.", "",
               "## Fixed pedigree with more genotyped animals", "",
               markdown(geno, [("n", "Pedigree size"), ("genotipados", "Genotyped animals"), ("operador", "Operator"),
                               ("tiempo_pca_min", "PCA (min)"), ("rss_proceso_mib", "RSS (MiB)"),
                               ("residuo_max_relativo", "Relative residual")]), "",
               "q=4, l=20, with one PCA seed. The observed decrease in runtime does not establish that adding genotyped animals generally speeds up the method: repeated measurements and control of the computing environment are needed. Changing which animals are genotyped also changes the operator; the residual is not error against an exact PCA reference.", "",
               "The transferred animal-ID files confirm nested subsets. The data generator checked that pedigrees and shared genotypes were unchanged; this has not been rechecked locally because the input .rds datasets were not downloaded.", "",
               "## Run coverage", "",
               markdown([{**r, "coverage_en": "complete" if r["cobertura"] == "completa" else "incomplete"} for r in inventory],
                        [("run", "Run"), ("pca_esperadas", "Expected PCA runs"), ("pca_ok", "Successful PCA runs"),
                         ("diagnosticos_ok", "Successful diagnostics"), ("coverage_en", "Coverage")]), "",
               "No CHC results for 100,000 individuals are present in the downloaded files. These CSVs alone do not establish whether those two configurations were submitted or what happened to them.", "",
               "The explicit timeout is retained in the failure table and is not treated as a completed PCA runtime. A scheduler termination may leave no failure row in results.csv; accounting records or logs are needed to establish the cause. Scheduler accounting shared in this conversation confirmed the 48-hour termination of large_20261029_122040_658611.", "",
               "## Units and interpretation", "",
               "- PCA time: internal computation time, including pedigree factorization and Hutch++, excluding data loading, saving and subsequent diagnostics.",
               "- Process time: duration of the entire PCA process.",
               "- RSS (MiB or GiB): peak resident memory of the entire process, not incremental memory used exclusively by PCA.",
               "- Relative residual: max_j ||K u_j - d_j u_j|| / ||K u_j||. It is not relative error against exact eigenvalues.",
               "- Variance proportions depend on the Hutch++ trace estimate. Accurate eigenvalues do not guarantee equally accurate variance proportions.",
               "- Source runs and code signatures are retained for auditing. The available information is insufficient to attribute runtime differences to a particular node-level cause.", "",
               "## Reproduce the report", "", "Run these commands from the repository root on the Mac:", "", "```sh",
               "python3 validation/summarise_eddie.py",
               "Rscript validation/export_eddie_tables.R reports/eddie_comparison en",
               "python3 validation/render_eddie_report.py --language en", "```", "",
               "The summary uses only the Python standard library. Excel and figure export requires the R packages openxlsx and ggplot2, which are available on this computer.", ""]
    if args.language in ("en", "both"):
        (data_dir / "comparison_en.md").write_text("\n".join(english))
    print(f"Generated {len(scale)} scalability rows, {len(exact)} accuracy summaries, "
          f"{len(geno)} genotyping rows; {len(failures)} explicit failed stage(s).")
    print(f"Report data: {data_dir}")


if __name__ == "__main__":
    main()
