#!/usr/bin/env python3
"""
scVIRGILS - Tkinter GUI for running a Snakemake QC and Filtering pipeline
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
    "Modeling": {"job_id": None, "running": False, "status_text": "Waiting for Modeling inputs..."}
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
# Enable/Disable Run button based on readiness
# --------------------------
def check_all_ready():
    global job_running
    if cellranger_saved.get() and metadata_saved.get() and sample_key_saved.get() and seq_batch_key_saved.get():
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
    if stage in ["QC", "Filtering", "Modeling"]:
        sbatch_script = "snakemake.sh"
    else:
        return

    
    button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else model_run
    progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else model_progress_bar
    label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else model_status_label

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
# Enable Modeling stage
# --------------------------

def enable_modeling_stage():
    patch_snakefile_for_modeling()
    model_status_label.config(text="Modeling stage enabled.", fg='blue')
    move_to_modeling.config(state='disabled')  # disable after click

# --------------------------
# Check SLURM job status
# --------------------------
def check_job_status(job_id, stage="QC"):
    global job_state

    button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else model_run
    progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else model_progress_bar
    label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else model_status_label

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
            move_to_modeling.config(state='normal')

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
                    f.write("        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'")
                    inside_all = True
                elif inside_all and (line.strip().startswith("input:") or line.strip().startswith("#") or line.startswith(" ")):
                    continue
                else:
                    f.write(line)
                    inside_all = False
    except Exception as e:
        QC_status_label.config(text=f"Error patching Snakefile: {e}", fg='red')

# --------------------------
# Patch Snakefile for Modeling
# --------------------------
def patch_snakefile_for_modeling():
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
                    f.write("        merged_rna_anndata = work_dir+'/atlas/04_annotated_anndata_rna.h5ad'\n")
                    inside_all = True
                elif inside_all and (line.strip().startswith("input:") or line.strip().startswith("#") or line.startswith(" ")):
                    continue
                else:
                    f.write(line)
                    inside_all = False
    except Exception as e:
        model_status_label.config(text=f"Error patching Snakefile: {e}", fg='red')

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

def load_gui_state():
    global job_state
    if not os.path.exists(STATE_FILE):
        return
    try:
        with open(STATE_FILE, "r") as f:
            state = json.load(f)
    except Exception as e:
        QC_status_label.config(text=f"Error loading GUI state: {e}", fg='red')
        return

    job_state = state.get("job_state", job_state)

    for stage, info in job_state.items():
        button = QC_run if stage == "QC" else filter_run if stage == "Filtering" else model_run
        progress = QC_progressbar if stage == "QC" else filter_progressbar if stage == "Filtering" else model_progress_bar
        label = QC_status_label if stage == "QC" else filter_status_label if stage == "Filtering" else model_status_label

        label.config(text=info.get("status_text", label.cget("text")))
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

    # Auto-enable modeling if final filtered merged anndata exists
    if os.path.exists("atlas/02_filtered_anndata_rna.h5ad"):
        filter_status_label.config(text="Filtering complete (detected atlas).", fg="green")
        move_to_modeling.config(state='normal')


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
# --------------------------
entries = [
    ("CELLRANGER Path", "data_dir", cellranger_saved),
    ("METADATA Path", "metadata_table", metadata_saved),
    ("Sample Key (e.g. sample_id)", "sample_key", sample_key_saved),
    ("Seq Batch Key (e.g. sequencing_round)", "seq_batch_key", seq_batch_key_saved)
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
    elif var_name == "seq_batch_key":
        seq_batch_entry = entry

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

move_to_modeling = ttk.Button(
    scrollable_frame,
    text="Move to Modeling",
    command=enable_modeling_stage,
    state='disabled'
)

move_to_modeling.grid(row=22, column=0, columnspan=2, pady=15, sticky='w')

# --------------------------
# Interface Text / Header (Modeling)
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Modeling', style='Header.Label')
header.grid(row=23, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# GUI widgets for Modeling (after defining scrollable_frame)
# --------------------------
model_run = ttk.Button(scrollable_frame, text='Run Modeling!', command=lambda: start_snakemake_job(stage="Modeling"), state='disabled')
model_run.grid(row=24, column=3, pady=10, padx=10, sticky='w')

model_progress_bar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
model_progress_bar.grid(row=24, column=0, columnspan=3, padx=10, sticky='ew')

model_status_label = tk.Label(scrollable_frame, text="Waiting to start Modeling...", fg="black")
model_status_label.grid(row=25, column=0, columnspan=4, pady=10, sticky='w')



# --------------------------
# Load any previous GUI state after UI creation
# --------------------------
load_gui_state()

# --------------------------
# Start the GUI event loop
# --------------------------
root.mainloop()
