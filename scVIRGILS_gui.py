#!/usr/bin/env python3
"""
scVIRGILS - Tkinter GUI for running a Snakemake QC, Filtering, and batch correction pipeline
"""

import tkinter as tk
from tkinter import ttk
import subprocess
import re
from PIL import Image, ImageTk
import webbrowser
import json
import os

# --------------------------
# Globals and persistent file
# --------------------------
current_job_id = None
job_running = False
STATE_FILE = "gui_state.json"
# Conversion I/O
CONVERT_IN = os.path.join("atlas", "03_modeled_anndata_rna.h5ad")
CONVERT_OUT = os.path.join("atlas", "03_modeled_anndata_rna.h5seurat")
conversion_proc = None  # subprocess handle

CONVERT_SCRIPT = "/scripts/convert_h5ad_to_seurat.sh"
CONVERT_POLL_MS = 1500

qc_flags = {
    "mito": False,
    "ribo": False,
    "gene": False,
    "doublet": False,
    "genes_by_counts": False
}

job_state = {
    "QC": {"job_id": None, "running": False, "status_text": "Waiting for QC inputs..."},
    "Filtering": {"job_id": None, "running": False, "status_text": "Waiting for Filtering inputs..."},
    "Batch Correction": {"job_id": None, "running": False, "status_text": "Waiting for Batch Correction inputs..."},
    "Conversion": {"job_id": None, "running": False, "status_text": "Waiting for Pipeline completion..."}
}

# --------------------------
# Basic Tk root + layout
# --------------------------
root = tk.Tk()
root.geometry('1200x700')
root.minsize(900, 600)
root.title('scVIRGILS')

container = tk.Frame(root)
container.pack(fill=tk.BOTH, expand=True)

canvas = tk.Canvas(container, borderwidth=0)
canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
canvas.configure(yscrollcommand=scrollbar.set)

scrollable_frame = tk.Frame(canvas)
scrollable_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

# ------------ scrolling handlers ------------
def _on_mousewheel(event):
    canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

def on_frame_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

def on_canvas_configure(event):
    canvas.itemconfig(scrollable_window, width=event.width)

scrollable_frame.bind("<Configure>", on_frame_configure)
canvas.bind("<Configure>", on_canvas_configure)
root.bind_all("<MouseWheel>", _on_mousewheel)

# --------------------------
# Input state variables
# --------------------------
cellranger_saved = tk.BooleanVar(value=False)
metadata_saved = tk.BooleanVar(value=False)
sample_key_saved = tk.BooleanVar(value=False)
seq_batch_key_saved = tk.BooleanVar(value=False)

mito_percent_thresh = tk.StringVar(value="")
ribo_percent_thresh = tk.StringVar(value="")
doublet_thresh = tk.StringVar(value="")
min_genes_per_cell = tk.StringVar(value="")

# --------------------------
# Styling
# --------------------------
style = ttk.Style()
style.configure('Header.Label', font=('Helvetica', 16, 'bold'))
style.map('TButton',
          foreground=[('disabled', '#808080')],
          background=[('disabled', '#d3d3d3')])

# --------------------------
# Helper: write value to Snakefile
# --------------------------
def fill(variable_name, entered_path, status_label=None, flag_var=None, numeric=False):
    target_file = "snakefile"

    # Validate numeric inputs
    if numeric:
        try:
            float(entered_path)  # allow ints and floats
        except ValueError:
            if status_label:
                status_label.config(text="Not a valid numeric entry", fg='red')
            check_filter_ready()
            return
        new_line = f'{variable_name} = {entered_path}\n'
    else:
        new_line = f'{variable_name} = "{entered_path}"\n'

    replaced_once = False

    try:
        with open(target_file, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        if status_label:
            status_label.config(text=f"Error: File '{target_file}' not found.", fg='red')
        return

    try:
        with open(target_file, 'w') as file:
            for line in lines:
                if not replaced_once and line.strip().startswith(f'{variable_name} ='):
                    file.write(new_line)
                    replaced_once = True
                else:
                    file.write(line)

            if not replaced_once:
                file.write(new_line)
    except Exception as e:
        if status_label:
            status_label.config(text=f"Error writing '{target_file}': {e}", fg='red')
        return

    if status_label:
        status_label.config(text=f"{variable_name} saved!", fg='green')
    if flag_var:
        try:
            flag_var.set(True)
        except Exception:
            pass

    check_all_ready()
    check_filter_ready()

# --------------------------
# Conversion Helpers 
# --------------------------
def update_conversion_button_state():
    """Enable/disable Convert! based on presence of input file, script, and running state."""
    exists_in = os.path.exists(CONVERT_IN)
    script_ok = os.path.exists(CONVERT_SCRIPT)
    running = job_state["Conversion"].get("running", False)
    enable = exists_in and script_ok and not running
    convert_run.config(state=('normal' if enable else 'disabled'))
    root.after(3000, update_conversion_button_state)

def start_conversion_job():
    """Run the bash conversion script to turn h5ad -> h5seurat."""
    global conversion_proc

    if not os.path.exists(CONVERT_IN):
        convert_status_label.config(text=f"Input not found: {CONVERT_IN}", fg="red")
        return

    script_path = "/scripts/convert_h5ad_to_seurat.sh"
    if not os.path.exists(script_path):
        convert_status_label.config(text=f"Script not found: {script_path}", fg="red")
        return
    if not os.access(script_path, os.X_OK):
        try:
            os.chmod(script_path, 0o755)
        except Exception:
            convert_status_label.config(text=f"Script not executable (chmod +x): {script_path}", fg="red")
            return

    try:
        convert_run.config(state='disabled')
        convert_progress_bar.start()
        convert_status_label.config(text="Conversion started…", fg="blue")
        job_state["Conversion"].update({"running": True, "status_text": "Conversion in progress..."})
        save_gui_state()

        # Use a login shell so 'module' is available inside the script if needed.
        # No args: the script handles paths itself.
        conversion_proc = subprocess.Popen(
            ["/bin/bash", "-lc", script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=os.path.abspath(os.getcwd())  # run from repo root; script is robust anyway
        )

        root.after(1500, check_conversion_status)
    except Exception as e:
        convert_status_label.config(text=f"Error starting conversion: {e}", fg="red")
        convert_run.config(state='normal')


def check_conversion_status():
    """Consider conversion complete when CONVERT_OUT exists; otherwise handle failures."""
    global conversion_proc

    # Success path: output file appeared
    if os.path.exists(CONVERT_OUT):
        convert_progress_bar.stop()
        convert_progress_bar['value'] = 100
        job_state["Conversion"].update({"running": False, "status_text": "Conversion complete!", "progress": 100})
        convert_status_label.config(text="Conversion complete!", fg="green")
        convert_run.config(state='normal')
        save_gui_state()
        return

    # If process ended but no output, show error
    if conversion_proc is not None:
        rc = conversion_proc.poll()
        if rc is not None:
            out, err = conversion_proc.communicate()
            job_state["Conversion"].update({"running": False, "status_text": f"Conversion failed (code {rc})"})
            convert_progress_bar.stop()
            convert_status_label.config(text=f"Conversion failed (exit {rc}). See stderr in console.", fg="red")
            if out: print("[convert stdout]\n", out)
            if err: print("[convert stderr]\n", err)
            convert_run.config(state='normal')
            save_gui_state()
            return

    # Still running; poll again
    root.after(CONVERT_POLL_MS, check_conversion_status)

# --------------------------
# Enable/Disable Run button based on readiness
# --------------------------

def check_all_ready():
    """QC run is enabled once CORE inputs are saved. Seq batch key is no longer required here."""
    global job_running
    if cellranger_saved.get() and metadata_saved.get() and sample_key_saved.get():
        if not job_running:
            QC_run.config(state='normal')
        else:
            QC_run.config(state='disabled')
    else:
        QC_run.config(state='disabled')


def check_filter_ready():
    global job_running
    all_filled = all([
        mito_percent_thresh.get(),
        ribo_percent_thresh.get(),
        doublet_thresh.get(),
        min_genes_per_cell.get()
    ])
    numeric_invalid = False
    for val in [mito_percent_thresh.get(), ribo_percent_thresh.get(), doublet_thresh.get(), min_genes_per_cell.get()]:
        try:
            float(val)
        except ValueError:
            numeric_invalid = True
            break
    if all_filled and not numeric_invalid and not job_running:
        filter_run.config(state='normal')
    else:
        filter_run.config(state='disabled')

# --------------------------
# Start a Snakemake job via sbatch
# --------------------------

def start_snakemake_job(stage="QC"):
    global job_state

    # Explicit sbatch script for each stage
    if stage in ["QC", "Filtering", "Batch Correction"]:
        sbatch_script = "snakemake.sh"
    else:
        return

    
    button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else batch_correction_run if stage == "Batch Correction" else convert_run
    progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else batch_correction_progress_bar if stage == "Batch Correction" else convert_progress_bar
    label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else batch_correction_status_label if stage == "Batch Correction" else convert_status_label

    try:
        result = subprocess.run(['sbatch', sbatch_script], capture_output=True, check=True, text=True)
        stdout = result.stdout.strip()
        match = re.search(r'(\d+)', stdout)
        if match:
            job_id = match.group(1)
            job_state[stage].update({
                "job_id": job_id,
                "running": True,
                "status_text": f"{stage} job {job_id} submitted.",
                "progress": None
            })
            button.config(state='disabled')
            progress.start()
            label.config(text=job_state[stage]["status_text"], fg='blue')
            save_gui_state()
            root.after(200, lambda: check_job_status(job_id, stage))
        else:
            label.config(text=f"No job ID found in sbatch output: {stdout}", fg='red')
    except subprocess.CalledProcessError as cpe:
        stderr = cpe.stderr.strip() if cpe.stderr else str(cpe)
        label.config(text=f"sbatch error: {stderr}", fg='red')
        button.config(state='normal')
    except Exception as e:
        label.config(text=f"Error starting job: {e}", fg='red')
        button.config(state='normal')

# --------------------------
# Enable filtering stage
# --------------------------

def enable_filtering_stage():
    for btn in save_buttons.values():
        btn.config(state='normal')
    patch_snakefile_for_filtering()
    QC_status_label.config(text="Filtering stage enabled.", fg='blue')
    move_to_filtering.config(state='disabled')  # disable after click


# --------------------------
# Enable Batch Correction stage
# --------------------------

def enable_batch_correction_stage():
    patch_snakefile_for_batch_correction()
    batch_correction_status_label.config(text="Batch correction stage enabled.", fg='blue')
    move_to_batch_correction.config(state='disabled')  # disable after click

    # Enable the Seq Batch Key save button now that we're entering this stage
    try:
        seq_batch_save_button.config(state='normal')
    except Exception:
        pass

    if not job_state["Batch Correction"]["running"]:
        batch_correction_run.config(state='normal')
        
# --------------------------
# Check SLURM job status
# --------------------------

def check_job_status(job_id, stage="QC"):
    global job_state

    button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else batch_correction_run if stage == "Batch Correction" else convert_run
    progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else batch_correction_progress_bar if stage == "Batch Correction" else convert_progress_bar
    label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else batch_correction_status_label if stage == "Batch Correction" else convert_status_label

    try:
        result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
        stdout = result.stdout
    except Exception as e:
        label.config(text=f"Error checking job status: {e}", fg='red')
        root.after(5000, lambda: check_job_status(job_id, stage))
        return

    if job_id in stdout:
        job_state[stage]["running"] = True
        job_state[stage]["status_text"] = f"{stage} job {job_id} is still running..."
        label.config(text=job_state[stage]["status_text"], fg='orange')
        root.after(5000, lambda: check_job_status(job_id, stage))
    else:
        job_state[stage]["running"] = False
        job_state[stage]["status_text"] = f"{stage} job {job_id} is complete!"
        job_state[stage]["progress"] = 100
        label.config(text=job_state[stage]["status_text"], fg='green')
        progress.stop()
        progress['value'] = 100
        button.config(state='normal')
        save_gui_state()
        if stage == "QC":
            next_process()
        check_filter_ready()
        if stage == "Filtering":
            move_to_batch_correction.config(state='normal')
            next_process()
        elif stage == "Batch Correction":
            convert_status_label.config(text="Batch correction complete.", fg="green")
            # Button will auto-enable via update_conversion_button_state() within 3s


# --------------------------
# Enable QC result buttons
# --------------------------

def next_process():
    open_mito_qc.config(state='normal')
    open_ribo_qc.config(state='normal')
    open_gene_qc.config(state='normal')
    open_doublet_qc.config(state='normal')
    open_genes_by_counts_qc.config(state='normal')

def mark_qc_viewed(flag_key):
    qc_flags[flag_key] = True
    if all(qc_flags.values()):
        move_to_filtering.config(state='normal')

# --------------------------
# Patch Snakefile for filtering
# --------------------------

def patch_snakefile_for_filtering():
    target_file = "snakefile"
    try:
        with open(target_file, "r") as f:
            lines = f.readlines()
        with open(target_file, "w") as f:
            inside_all = False
            for line in lines:
                if line.strip().startswith("rule all"):
                    f.write("rule all:\n")
                    f.write("    input:\n")
                    f.write("        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'\n")

                    inside_all = True
                elif inside_all and (line.strip().startswith("input:") or line.strip().startswith("#") or line.startswith(" ")):
                    continue
                else:
                    f.write(line)
                    inside_all = False
    except Exception as e:
        QC_status_label.config(text=f"Error patching Snakefile: {e}", fg='red')

# --------------------------
# Patch Snakefile for Batch Correction
# --------------------------

def patch_snakefile_for_batch_correction():
    target_file = "snakefile"
    try:
        with open(target_file, "r") as f:
            lines = f.readlines()

        with open(target_file, "w") as f:
            inside_all = False
            for line in lines:
                stripped = line.strip()

                # Detect the start of rule all
                if stripped.startswith("rule all"):
                    inside_all = True
                    f.write(line)
                    continue

                # Detect input line inside rule all
                if inside_all and "merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'" in stripped:
                    f.write("        merged_rna_anndata = work_dir+'/atlas/03_modeled_anndata_rna.h5ad'\n")
                    continue

                # Stop tracking once we hit a """ or something not part of the input
                if inside_all and stripped.startswith('"""'):
                    f.write(line)
                    inside_all = False
                    continue

                # Write everything else as-is
                f.write(line)

    except Exception as e:
        batch_correction_status_label.config(text=f"Error patching Snakefile: {e}", fg='red')


# --------------------------
# GUI state persistence
# --------------------------

def save_gui_state(stage=None):
    state = {
        "data_dir": data_dir_entry.get() if 'data_dir_entry' in globals() else "",
        "metadata_table": metadata_dir_entry.get() if 'metadata_dir_entry' in globals() else "",
        "sample_key": sample_key_entry.get() if 'sample_key_entry' in globals() else "",
        "seq_batch_key": seq_batch_entry.get() if 'seq_batch_entry' in globals() else "",
        "cellranger_saved": cellranger_saved.get(),
        "metadata_saved": metadata_saved.get(),
        "sample_key_saved": sample_key_saved.get(),
        "seq_batch_key_saved": seq_batch_key_saved.get(),
        "mito_percent_thresh": mito_percent_thresh.get(),
        "ribo_percent_thresh": ribo_percent_thresh.get(),
        "doublet_thresh": doublet_thresh.get(),
        "min_genes_per_cell": min_genes_per_cell.get(),
        "job_state": job_state
    }
    try:
        with open(STATE_FILE, "w") as f:
            json.dump(state, f)
    except Exception as e:
        QC_status_label.config(text=f"Error saving GUI state: {e}", fg='red')


for stage, info in job_state.items():
    button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else batch_correction_run if stage == "Batch Correction" else convert_run
    progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else batch_correction_progress_bar if stage == "Batch Correction" else convert_progress_bar
    label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else batch_correction_status_label if stage == "Batch Correction" else convert_status_label

    # Always restore the last status text
    label.config(text=info.get("status_text", label.cget("text")))

    if stage == "Conversion":
        # Do NOT auto-start the Convert progress bar.
        # If a previous session marked it running, keep button disabled and just poll for completion.
        if info.get("running", False):
            button.config(state='disabled')
            # Resume polling BUT leave the bar idle until the user actually clicks Convert! in this session.
            root.after(CONVERT_POLL_MS, check_conversion_status)
        elif info.get("progress") is not None:
            progress.stop()
            progress['value'] = info["progress"]
            button.config(state='normal')
        continue  # skip the generic auto-start logic below

    # Generic resume behavior for QC / Filtering / Batch Correction
    if info.get("running", False):
        progress.start()
        button.config(state='disabled')
        if info.get("job_id"):
            root.after(200, lambda j=info.get("job_id"), s=stage: check_job_status(j, s))
    elif info.get("progress") is not None:
        progress.stop()
        progress['value'] = info["progress"]
        button.config(state='normal')

    # Auto-enable filtering if QC result file exists
    if os.path.exists("figures/QC_mito_pct.png"):
        QC_status_label.config(text="QC complete (detected QC results).", fg="green")
        next_process()  # enable QC view buttons
        move_to_filtering.config(state='normal')

    # Auto-enable batch_correction if final filtered merged anndata exists
    if os.path.exists("atlas/02_filtered_anndata_rna.h5ad"):
        filter_status_label.config(text="Filtering complete (detected atlas).", fg="green")
        move_to_batch_correction.config(state='normal')

    # Auto-enable convert if batch corrected data exists
    # Auto-enable convert if modeled anndata exists
    if os.path.exists(CONVERT_IN):
        convert_status_label.config(text="Batch correction complete (detected atlas).", fg="green")
        # Only enable if not currently running
        if not job_state["Conversion"].get("running", False):
            convert_run.config(state='normal')

    # If output already exists, show as complete
    if os.path.exists(CONVERT_OUT):
        convert_status_label.config(text="Conversion complete (detected h5seurat).", fg="green")
        convert_run.config(state='normal')


# --------------------------
# Handle window close
# --------------------------
root.protocol("WM_DELETE_WINDOW", lambda: (save_gui_state(), root.destroy()))


# --------------------------
# Interface Text / Header
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Quality Control', style='Header.Label')
header.grid(row=0, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Optional Image (logo)
# --------------------------
try:
    img = Image.open("images/VIRGIL.png")
    img = img.resize((150, 150))
    photo = ImageTk.PhotoImage(img)
    logo_label = tk.Label(scrollable_frame, image=photo)
    logo_label.image = photo
    logo_label.grid(row=0, column=4, padx=10, sticky='e')
except Exception:
    pass

# --------------------------
# Primary Input Entries (first block)
# NOTE: Seq Batch Key is REMOVED from this block and added near Batch Correction
# --------------------------
entries = [
    ("CELLRANGER Path", "data_dir", cellranger_saved),
    ("METADATA Path", "metadata_table", metadata_saved),
    ("Sample Key (e.g. sample_id)", "sample_key", sample_key_saved)
]

for i, (label_text, var_name, flag_var) in enumerate(entries):
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=i+1, column=0, sticky='e', padx=10, pady=5)
    entry = ttk.Entry(scrollable_frame, width=50)
    entry.grid(row=i+1, column=1, sticky='w')
    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=i+1, column=3, sticky='w')
    btn = ttk.Button(scrollable_frame, text='Save', command=lambda v=var_name, e=entry, s=status, f=flag_var: fill(v, e.get(), s, f))
    btn.grid(row=i+1, column=2, sticky='w', padx=5)

    if var_name == "data_dir":
        data_dir_entry = entry
    elif var_name == "metadata_table":
        metadata_dir_entry = entry
    elif var_name == "sample_key":
        sample_key_entry = entry

# --------------------------
# Run + Progress + Status (QC)
# --------------------------
QC_run = ttk.Button(scrollable_frame, text='Run QC!', command=start_snakemake_job, state='disabled')
QC_run.grid(row=6, column=3, pady=10, padx=10, sticky='w')

QC_progressbar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
QC_progressbar.grid(row=6, column=0, columnspan=3, padx=10, sticky='ew')

QC_status_label = tk.Label(scrollable_frame, text="Waiting for QC inputs...", fg="black")
QC_status_label.grid(row=7, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Interface Text / Header (Filtering)
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Filtering', style='Header.Label')
header.grid(row=9, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# QC View buttons
# --------------------------
qc_key_map = {
    "% Mitochondria": "mito",
    "% Ribosomal": "ribo",
    "Gene Counts": "gene",
    "Doublet Score": "doublet",
    "Genes by Counts": "genes_by_counts"
}

view_buttons = [
    ("% Mitochondria", 'figures/QC_mito_pct.png'),
    ("% Ribosomal", 'figures/QC_ribo_pct.png'),
    ("Gene Counts", 'figures/QC_gene_counts.png'),
    ("Doublet Score", 'figures/QC_doublet.png'),
    ("Genes by Counts", 'figures/QC_genes_by_counts.png')
]

for i, (text, path) in enumerate(view_buttons):
    btn = ttk.Button(
        scrollable_frame,
        text=f"View {text} QC",
        command=lambda p=path, k=qc_key_map[text]: (webbrowser.open(f"file://{os.path.abspath(p)}"), mark_qc_viewed(k)),
        state='disabled'
    )
    btn.grid(row=10+i, column=0, columnspan=2, padx=10, pady=5, sticky='w')

    if i == 0:
        open_mito_qc = btn
    elif i == 1:
        open_ribo_qc = btn
    elif i == 2:
        open_gene_qc = btn
    elif i == 3:
        open_doublet_qc = btn
    elif i == 4:
        open_genes_by_counts_qc = btn

move_to_filtering = ttk.Button(
    scrollable_frame,
    text="Move to Filtering",
    command=enable_filtering_stage,
    state='disabled'
)
move_to_filtering.grid(row=15, column=0, columnspan=2, pady=15, sticky='w')

# --------------------------
# Threshold entries (Filtering)
# --------------------------
base_row = 16

threshold_entries = [
    ("Mitochondria % threshold (e.g. 20)", "mito_percent_thresh", mito_percent_thresh),
    ("Ribosomal % threshold (e.g. 20)", "ribo_percent_thresh", ribo_percent_thresh),
    ("Doublet threshold (e.g. 0.15)", "doublet_thresh", doublet_thresh),
    ("Minimum genes per cell (e.g. 200)", "min_genes_per_cell", min_genes_per_cell)
]

entry_widgets = {}
save_buttons = {}

for i, (label_text, var_name, tk_var) in enumerate(threshold_entries):
    row_idx = base_row + i
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=row_idx, column=0, sticky='e', padx=10, pady=5)

    entry = ttk.Entry(scrollable_frame, width=50, textvariable=tk_var)
    entry.grid(row=row_idx, column=1, sticky='w')

    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=row_idx, column=3, sticky='w')

    btn = ttk.Button(
        scrollable_frame,
        text='Save',
        command=lambda v=var_name, e=entry, s=status: (fill(v, e.get(), s, None, numeric=True), check_filter_ready()),
        state='disabled'
    )
    btn.grid(row=row_idx, column=2, sticky='w', padx=5)

    entry_widgets[var_name] = entry
    save_buttons[var_name] = btn

# --------------------------
# Run + Progress + Status (Filtering)
# --------------------------
filter_run = ttk.Button(scrollable_frame, text='Run Filtering!', command=lambda: start_snakemake_job(stage="Filtering"), state='disabled')
filter_run.grid(row=20, column=3, pady=10, padx=10, sticky='w')

filter_progressbar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
filter_progressbar.grid(row=20, column=0, columnspan=3, padx=10, sticky='ew')

filter_status_label = tk.Label(scrollable_frame, text="Waiting for Filtering inputs...", fg="black")
filter_status_label.grid(row=21, column=0, columnspan=4, pady=10, sticky='w')

move_to_batch_correction = ttk.Button(
    scrollable_frame,
    text="Move to Batch correction",
    command=enable_batch_correction_stage,
    state='disabled'
)

move_to_batch_correction.grid(row=22, column=0, columnspan=2, pady=15, sticky='w')

# --------------------------
# Seq Batch Key input (moved here, just before Batch Correction)
# Save button is DISABLED until "Move to Batch correction" is clicked
# --------------------------
seq_batch_row = 22  # place directly after the move_to_batch_correction button

seq_batch_label = ttk.Label(scrollable_frame, text="Seq Batch Key (e.g. sequencing_round)", wraplength=400)
seq_batch_label.grid(row=seq_batch_row+1, column=0, sticky='e', padx=10, pady=5)

seq_batch_entry = ttk.Entry(scrollable_frame, width=50)
seq_batch_entry.grid(row=seq_batch_row+1, column=1, sticky='w')

seq_batch_status = tk.Label(scrollable_frame, text="", anchor='w')
seq_batch_status.grid(row=seq_batch_row+1, column=3, sticky='w')

seq_batch_save_button = ttk.Button(
    scrollable_frame,
    text='Save',
    command=lambda: fill("seq_batch_key", seq_batch_entry.get(), seq_batch_status, seq_batch_key_saved),
    state='disabled'  # <- remains disabled until user clicks Move to Batch correction
)
seq_batch_save_button.grid(row=seq_batch_row+1, column=2, sticky='w', padx=5)

# --------------------------
# Interface Text / Header (Batch correction)
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Batch Correction', style='Header.Label')
header.grid(row=23+1, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# GUI widgets for Batch Correction (after defining scrollable_frame)
# --------------------------
batch_correction_run = ttk.Button(scrollable_frame, text='Run Batch Correction!', command=lambda: start_snakemake_job(stage="Batch Correction"), state='disabled')
# placed after header row
batch_correction_run.grid(row=24+1, column=3, pady=10, padx=10, sticky='w')

batch_correction_progress_bar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
batch_correction_progress_bar.grid(row=24+1, column=0, columnspan=3, padx=10, sticky='ew')

batch_correction_status_label = tk.Label(scrollable_frame, text="Waiting to start Batch Correction...", fg="black")
batch_correction_status_label.grid(row=25+1, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Interface Text / Header (Convert to seurat)
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Convert h5ad to seurat', style='Header.Label')
header.grid(row=26+1, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# GUI widgets for converting to seurat object (after defining scrollable_frame)
# --------------------------
convert_run = ttk.Button(scrollable_frame, text='Convert!', command=start_conversion_job, state='disabled')
# placed after header row
convert_run.grid(row=27+1, column=3, pady=10, padx=10, sticky='w')

convert_progress_bar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
convert_progress_bar.grid(row=27+1, column=0, columnspan=3, padx=10, sticky='ew')

convert_status_label = tk.Label(scrollable_frame, text="Waiting to start Conversion...", fg="black")
convert_status_label.grid(row=28+1, column=0, columnspan=4, pady=10, sticky='w')


update_conversion_button_state()

# --------------------------
# Load any previous GUI state after UI creation
# --------------------------
load_gui_state()

# --------------------------
# Start the GUI event loop
# --------------------------
root.mainloop()
