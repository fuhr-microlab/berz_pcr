import streamlit as st
from pathlib import Path
import subprocess
import sys
import pandas as pd
import shutil

# Run with:  streamlit run app.py

import base64



def set_bg(image_path):
    with open(image_path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode()
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/png;base64,{encoded}");
            background-size: cover;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

set_bg("background3.png")


st.set_page_config(page_title="Berz PCR Pipeline", layout="wide")
st.title("Berz PCR Pipeline: Shear → In-silico PCR → Thermocycle → Serotype Report")

ROOT = Path(".").resolve()

def run_cmd(cmd, cwd=None):
    """Run a command and stream output into the UI."""
    st.code(" ".join(cmd))
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    output_lines = []
    for line in p.stdout:
        output_lines.append(line.rstrip("\n"))
        st.text(line.rstrip("\n"))
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {p.returncode}")
    return "\n".join(output_lines)

def reset_dir(path: Path, label: str):
    if path.exists():
        st.warning(f"{label} exists — deleting: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


# -----------------------
# INPUTS
# -----------------------
st.subheader("1) Inputs")

col1, col2 = st.columns(2)

with col1:
    resources_dir = st.text_input("Input cps loci folder (FASTA)", value="resources")
    primers_bed = st.text_input("Primer BED file", value="primal_resources_5/primer.bed")
    refs_fasta = st.text_input("cps reference FASTA", value="cps_refs.fasta")

with col2:
    holocene_dir = st.text_input("Sheared reads output folder", value="holocene")

    # NEW: split raw products (cycle 0) vs cycled products (final)
    products_raw_dir = st.text_input("PCR products (raw) folder", value="products")
    products_dir = st.text_input("PCR products (cycled) folder", value="products_cycled")  # NEW default

    report_script = st.text_input("Report script", value="serotype_from_products.py")

# NEW: amplification script path
st.subheader("1b) Thermocycling script (optional)")
amplification_script = st.text_input("Amplification script", value="amplification.py")

st.subheader("2) Shearing settings (main.py)")
c1, c2, c3, c4 = st.columns(4)
with c1:
    chunk_size = st.number_input("Chunk size (bp)", min_value=50, max_value=5000, value=1000, step=50)
with c2:
    target_reads = st.number_input("Target reads", min_value=1, max_value=100000, value=2000, step=50)
with c3:
    seed = st.text_input("Random seed (blank = random each run)", value="")
with c4:
    avoid_dups = st.checkbox("Avoid duplicate windows", value=False)

st.subheader("2b) PCR mismatch settings")

max_mismatches = st.number_input(
    "Max mismatches per primer",
    min_value=0, max_value=10, value=2, step=1,
    key="max_mismatches"
)

three_prime_lock = st.number_input(
    "3′ end must match exactly (bases)",
    min_value=0, max_value=30, value=8, step=1,
    key="three_prime_lock"
)


# NEW: cycling settings
st.subheader("2b) Thermocycling settings (amplification.py)")
cc1, cc2, cc3 = st.columns(3)
with cc1:
    cycles = st.number_input("PCR cycles", min_value=1, max_value=60, value=35, step=1)
with cc2:
    n_out = st.number_input("Cycled output files (N)", min_value=1, max_value=200000, value=10000, step=100)
with cc3:
    amp_seed_text = st.text_input("Cycling RNG seed (leave blank for random)", value="")


st.subheader("3) Run controls")
run_shear = st.checkbox("Run shear (main.py)", value=True)
run_pcr = st.checkbox("Run in-silico PCR (insilico.py first round deterministic)", value=True)

# NEW
run_cycle = st.checkbox("Run thermocycling (amplification.py stochastic amplification)", value=True)

run_report = st.checkbox("Run report (serotype_from_products.py)", value=True)

#stoichastic
weighted_templates = st.checkbox(
    "Stochastic Gene Fragmentation (Requires very large target reads count if metagenome is large)",
    value=True
)


st.divider()

# -----------------------
# RUN
# -----------------------
if st.button("▶ Run Pipeline", type="primary"):
    try:
        # Validate paths (basic)
        resources_path = ROOT / resources_dir
        primers_path = ROOT / primers_bed
        refs_path = ROOT / refs_fasta

        if run_shear and not resources_path.exists():
            st.error(f"Input folder not found: {resources_path}")
            st.stop()
        if run_pcr and not primers_path.exists():
            st.error(f"Primer BED not found: {primers_path}")
            st.stop()
        if run_report and not refs_path.exists():
            st.error(f"Reference FASTA not found: {refs_path}")
            st.stop()

        st.success("Starting…")

        # 1) SHEAR
        if run_shear:
            st.subheader("Shear step (main.py)")

            holocene_path = ROOT / holocene_dir
            reset_dir(holocene_path, "Sheared reads output folder")

            cmd = [
                sys.executable, "main.py",
                "--input", str(resources_path),
                "--output", str(holocene_path),
                "--chunk", str(int(chunk_size)),
                "--target", str(int(target_reads)),
                "--seed", seed.strip() if seed.strip() else "none",
                "--avoid-dups", "1" if avoid_dups else "0",
                "--weighted", "1" if weighted_templates else "0",
            ]
            run_cmd(cmd, cwd=str(ROOT))

        # 2) PCR (raw products)
        if run_pcr:
            st.subheader("In-silico PCR (insilico.py)")

            products_raw_path = ROOT / products_raw_dir
            reset_dir(products_raw_path, "PCR products (raw) folder")

            cmd = [
                sys.executable, "insilico.py",
                "--templates", str(ROOT / holocene_dir),
                "--bed", str(primers_path),
                "--out", str(products_raw_path),

                "--max-mismatches", str(int(max_mismatches)),
                "--three-prime-lock", str(int(three_prime_lock)),
            ]

            run_cmd(cmd, cwd=str(ROOT))

            

        # 2b) Thermocycling (cycled products)
        if run_cycle:
            st.subheader("Thermocycling (amplification.py)")

            products_raw_path = ROOT / products_raw_dir
            
            report_path = products_raw_path / "amplification_report.tsv"
            if not report_path.exists():
                st.error(f"Missing report for cycling: {report_path}. Did insilico.py run?")
                st.stop()

            products_cycled_path = ROOT / products_dir
            reset_dir(products_cycled_path, "PCR products (cycled) folder")

            seed_args = []
            if amp_seed_text.strip():
                seed_args = ["--seed", amp_seed_text.strip()]

            cmd = [
                sys.executable, amplification_script,
                "--report", str(report_path),
                "--products-dir", str(products_raw_path),
                "--out-dir", str(products_cycled_path),
                "--cycles", str(int(cycles)),
                "--n-out", str(int(n_out)),
            ] + seed_args
            run_cmd(cmd, cwd=str(ROOT))

            png = products_cycled_path / "ct_curves.png"
            if png.exists():
                st.image(str(png), caption="Synthetic Ct curves (log fluorescence)")


        # 3) REPORT
        if run_report:
            st.subheader("Serotype report")

            # IMPORTANT: make report use cycled products by default
            # Option 1 (recommended): update serotype_from_products.py to accept --products
            # We'll call it that way if possible; if not, it will just run as-is.
            products_used = ROOT / (products_dir if run_cycle else products_raw_dir)

            cmd = [sys.executable, report_script]

            # If your report script supports CLI args, this is ideal:
            # e.g., serotype_from_products.py --products products_cycled --refs cps_refs.fasta
            cmd += ["--products", str(products_used)]  # safe if script accepts it
            cmd += ["--refs", str(refs_path)]          # safe if script accepts it

            try:
                run_cmd(cmd, cwd=str(ROOT))
            except RuntimeError:
                # Fallback: run without args (old behavior)
                st.warning("Report script did not accept --products/--refs; running without args.")
                run_cmd([sys.executable, report_script], cwd=str(ROOT))

            ranked = products_used / "serotype_ranked.tsv"
            if ranked.exists():
                df = pd.read_csv(ranked, sep="\t")
                st.subheader("Top serotype candidates")
                st.dataframe(df.head(20), width='stretch')
            else:
                st.warning("serotype_ranked.tsv not found (report may have failed or output path differs or no amplification occured with given sequences).")

        st.success("Pipeline finished.")

    except Exception as e:
        st.error(str(e))
