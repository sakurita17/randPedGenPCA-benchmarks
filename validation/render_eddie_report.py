"""Render the comparison Markdown as an offline, readable HTML report (stdlib)."""
import argparse
import base64
import html
from pathlib import Path
import re
import unicodedata


def inline(text):
    text = html.escape(text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    return text


def slug(text):
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", "-", text).strip("-")


def render(folder, language="en"):
    english = language == "en"
    data_dir = folder / "data" if (folder / "data").is_dir() else folder
    figures_dir = folder / "figures" if (folder / "figures").is_dir() else folder
    lines = (data_dir / ("comparison_en.md" if english else "comparacion.md")).read_text().splitlines()
    plots = {"Escalabilidad": "escalabilidad_tiempo.png",
             "Precisión frente a referencia exacta": "precision_eigenvalores.png",
             "Pedigrí fijo y más genotipados": "genotipados_tiempo.png"}
    technical = {"Cobertura de ejecuciones", "Unidades y límites", "Regenerar en el Mac"}
    if english:
        plots = {"Scalability": "escalabilidad_tiempo_en.png",
                 "Accuracy against an exact reference": "precision_eigenvalores_en.png",
                 "Fixed pedigree with more genotyped animals": "genotipados_tiempo_en.png"}
        technical = {"Run coverage", "Units and interpretation", "Reproduce the report"}
    parts = []
    i = 0
    section = None
    while i < len(lines):
        line = lines[i]
        if line.startswith("# "):
            i += 1
            continue
        if line.startswith("## "):
            if section:
                parts.append("</details>" if section == "details" else "</section>")
            title = line[3:]
            if title in technical:
                parts.append(f'<details id="{slug(title)}"><summary>{inline(title)}</summary>')
                section = "details"
            else:
                parts.append(f'<section id="{slug(title)}"><h2>{inline(title)}</h2>')
                section = "section"
            if title in plots and (figures_dir / plots[title]).exists():
                data = base64.b64encode((figures_dir / plots[title]).read_bytes()).decode()
                parts.append(f'<figure><img src="data:image/png;base64,{data}" '
                             f'alt="Gráfico de {html.escape(title)}"></figure>')
            i += 1
            continue
        if line.startswith("```"):
            block = []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                block.append(lines[i])
                i += 1
            parts.append("<pre><code>" + html.escape("\n".join(block)) + "</code></pre>")
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = [cell.strip() for cell in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", cell) for cell in cells):
                    rows.append(cells)
                i += 1
            parts.append('<div class="table-scroll" tabindex="0" role="region" aria-label="Tabla comparativa"><table>')
            parts.append("<thead><tr>" + "".join("<th scope=\"col\">" + inline(c) + "</th>" for c in rows[0]) + "</tr></thead><tbody>")
            for row in rows[1:]:
                parts.append("<tr>" + "".join("<td>" + inline(c) + "</td>" for c in row) + "</tr>")
            parts.append("</tbody></table></div>")
            continue
        if line.startswith("- "):
            parts.append("<ul>")
            while i < len(lines) and lines[i].startswith("- "):
                parts.append("<li>" + inline(lines[i][2:]) + "</li>")
                i += 1
            parts.append("</ul>")
            continue
        if line.strip():
            parts.append("<p>" + inline(line) + "</p>")
        i += 1
    if section:
        parts.append("</details>" if section == "details" else "</section>")
    document = '''<!doctype html>
<html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Resumen de los experimentos en Eddie</title>
<style>
:root{color-scheme:light;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;color:#17293b;background:#f2f5f8}
*{box-sizing:border-box}body{margin:0;line-height:1.65}main{max-width:1200px;margin:auto;padding:38px 24px 70px}
h1{font-size:clamp(28px,4vw,42px);line-height:1.15;margin:8px 0 18px}h2{font-size:26px;margin:0 0 16px}
.eyebrow{color:#425d73;letter-spacing:.1em;font-size:13px;font-weight:700;text-transform:uppercase}
.lead{font-size:19px;max-width:850px}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:24px 0}
.card,section,details{background:white;border:1px solid #d9e2ea;border-radius:12px;padding:24px}
.card strong{display:block;color:#145a77;font-size:22px}.card p{margin:6px 0 0}
nav{display:flex;gap:12px;flex-wrap:wrap;margin:26px 0}a{color:#145a77;text-underline-offset:3px}
nav a,.button{display:inline-block;padding:9px 15px;border:1px solid #b8cbd8;border-radius:8px;background:#fff}
section{margin:24px 0;scroll-margin-top:20px}p{margin:14px 0}li{margin:8px 0}
.table-scroll{overflow-x:auto;margin:20px 0;border:1px solid #d9e2ea;border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:15px;font-variant-numeric:tabular-nums}
th{background:#e7f0f5;text-align:left;white-space:normal;min-width:85px}th,td{padding:12px 14px;border-bottom:1px solid #e1e7ed}
td{white-space:nowrap}tbody tr:nth-child(even){background:#f7f9fb}tbody tr:hover{background:#ecf4f7}
figure{margin:16px 0}img{display:block;width:100%;height:auto}details{margin:14px 0}summary{cursor:pointer;font-size:19px;font-weight:650}
code{font-size:.9em;background:#eef2f6;padding:2px 5px;border-radius:4px;overflow-wrap:anywhere}pre{padding:16px;background:#eef2f6;overflow:auto}
.note{padding:16px 20px;border-left:4px solid #c18b30;background:#fff8e8;border-radius:5px}footer{color:#526779;margin-top:28px;font-size:14px}
@media(max-width:750px){main{padding:24px 14px}.cards{grid-template-columns:1fr}section,details{padding:18px}th,td{padding:10px}}
@media print{body{background:white}main{max-width:none;padding:0}.cards{break-inside:avoid}section{break-before:page}nav,.button{display:none}details{display:block}table{font-size:11px}figure{max-width:800px}section,details{border:0;padding:0}.table-scroll{overflow:visible}}
</style></head><body><main>
<header><div class="eyebrow">Resultados de randPedGenPCA · Eddie</div>
<h1>Resumen de los experimentos</h1>
<p class="lead">Empieza por los tres bloques principales: tamaño del pedigrí, precisión y número de genotipados. Cada bloque incluye un gráfico, su tabla y las limitaciones de la comparación.</p>
<div class="cards">
<article class="card"><strong>Hasta 100.000</strong><p>Individuos en el pedigrí. H tiene resultados completos en ese tamaño; faltan resultados descargados de CHC.</p></article>
<article class="card"><strong>Referencia exacta</strong><p>Comparación en 200, 500 y 1.000 individuos, con dos ajustes y tres semillas de PCA.</p></article>
<article class="card"><strong>Pedigrí fijo: 20.000</strong><p>Se comparan 500, 1.000, 2.000 y 4.000 genotipados de la misma población.</p></article>
</div>
<p><a class="button" href="cuadros_comparativos.xlsx">Abrir tablas en Excel</a></p>
<nav aria-label="Secciones del informe"><a href="#escalabilidad">1. Tamaño y coste</a><a href="#precision-frente-a-referencia-exacta">2. Precisión</a><a href="#pedigri-fijo-y-mas-genotipados">3. Más genotipados</a></nav>
<aside class="note"><strong>Cómo leer las tablas.</strong> H es el operador sin centrar y CHC es el centrado. q indica las iteraciones de potencia y l el número de vectores del sketch. PCA (h), (min) o (s) es tiempo interno de cálculo. RSS es la memoria máxima del proceso completo. Un residuo pequeño no equivale directamente a un error pequeño frente a la PCA exacta.</aside>
</header>
''' + "\n".join(parts) + '''
<footer>Informe local a partir de los CSV descargados de Eddie. Los gráficos están incluidos en este archivo: se puede leer sin conexión. Para actualizarlo, regenerar las tablas, los gráficos y después ejecutar validation/render_eddie_report.py.</footer>
</main></body></html>'''
    if english:
        translations = {
            'lang="es"': 'lang="en"',
            'Resumen de los experimentos en Eddie': 'Eddie experiment summary',
            'Resultados de randPedGenPCA · Eddie': 'randPedGenPCA results · Eddie',
            'Resumen de los experimentos': 'Experiment summary',
            'Empieza por los tres bloques principales: tamaño del pedigrí, precisión y número de genotipados. Cada bloque incluye un gráfico, su tabla y las limitaciones de la comparación.':
                'Start with the three main sections: pedigree size, accuracy and the number of genotyped animals. Each section includes a figure, a comparison table and guidance on interpretation.',
            'Hasta 100.000': 'Up to 100,000',
            'Individuos en el pedigrí. H tiene resultados completos en ese tamaño; faltan resultados descargados de CHC.':
                'Individuals in the pedigree. H has completed results at this size; CHC results are absent from the download.',
            'Referencia exacta': 'Exact reference',
            'Comparación en 200, 500 y 1.000 individuos, con dos ajustes y tres semillas de PCA.':
                'Comparisons at 200, 500 and 1,000 individuals, with two settings and three PCA seeds.',
            'Pedigrí fijo: 20.000': 'Fixed pedigree: 20,000',
            'Se comparan 500, 1.000, 2.000 y 4.000 genotipados de la misma población.':
                'Comparing 500, 1,000, 2,000 and 4,000 genotyped animals from the same population.',
            'cuadros_comparativos.xlsx': 'comparison_tables_en.xlsx',
            'Abrir tablas en Excel': 'Open comparison tables in Excel',
            'Secciones del informe': 'Report sections',
            '#escalabilidad': '#scalability',
            '#precision-frente-a-referencia-exacta': '#accuracy-against-an-exact-reference',
            '#pedigri-fijo-y-mas-genotipados': '#fixed-pedigree-with-more-genotyped-animals',
            '1. Tamaño y coste': '1. Size and computational cost',
            '2. Precisión': '2. Accuracy',
            '3. Más genotipados': '3. More genotyped animals',
            'Cómo leer las tablas.': 'How to read the tables.',
            'H es el operador sin centrar y CHC es el centrado. q indica las iteraciones de potencia y l el número de vectores del sketch. PCA (h), (min) o (s) es tiempo interno de cálculo. RSS es la memoria máxima del proceso completo. Un residuo pequeño no equivale directamente a un error pequeño frente a la PCA exacta.':
                'H is the uncentered operator and CHC is its centered version. q is the number of power iterations and l is the number of sketch vectors. PCA (h), (min) or (s) denotes internal computation time. RSS is peak resident memory of the entire process. A small residual does not directly imply a small error against exact PCA.',
            'Gráfico de ': 'Figure: ',
            'Tabla comparativa': 'Comparison table',
            'Informe local a partir de los CSV descargados de Eddie. Los gráficos están incluidos en este archivo: se puede leer sin conexión. Para actualizarlo, regenerar las tablas, los gráficos y después ejecutar validation/render_eddie_report.py.':
                'Local report generated from CSV files downloaded from Eddie. Figures are embedded, so the report can be read offline. To update it, regenerate the tables and English figures, then run validation/render_eddie_report.py --language en.'
        }
        for old, new in translations.items():
            document = document.replace(old, new)
    target = folder / ("summary_en.html" if english else "resumen.html")
    target.write_text(document)
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=Path("reports/eddie_comparison"))
    parser.add_argument("--language", choices=("es", "en"), default="en")
    args = parser.parse_args()
    print(render(args.folder, args.language))
