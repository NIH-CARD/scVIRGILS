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

# Track which QC buttons have been clicked
qc_flags = {
    "mito": False,
    "ribo": False,
    "gene": False,
    "doublet": False,
    "genes_by_counts": False
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
    if numeric:
        new_line = f'{variable_name} = {entered_path}\n'
    else:
        new_line = f'{variable_name} = "{entered_path}"\n'

    updated = False

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
                if line.strip().startswith(f'{variable_name} ='):
                    file.write(new_line)
                    updated = True
                else:
                    file.write(line)
            if not updated:
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
    """
    Enable the 'Run Filtering' button only when all threshold entries are saved and non-empty,
    and no filtering job is currently running.
    """
    global job_running
    all_filled = all([
        mito_percent_thresh.get(),
        ribo_percent_thresh.get(),
        doublet_thresh.get(),
        min_genes_per_cell.get()
    ])
    if all_filled and not job_running:
        filter_run.config(state='normal')
    else:
        filter_run.config(state='disabled')

# --------------------------
# Start a Snakemake job via sbatch
# --------------------------
def start_snakemake_job(stage="QC"):
    """
    stage: "QC" or "Filtering"
    Submits the appropriate sbatch script and tracks the job.
    """
    global current_job_id, job_running

    sbatch_script = "snakemake.sh" if stage == "QC" else "snakemake_filtering.sh"
    button = QC_run if stage == "QC" else filter_run

    try:
        result = subprocess.run(['sbatch', sbatch_script], capture_output=True, check=True, text=True)
        stdout = result.stdout.strip()
        match = re.search(r'(\d+)', stdout)
        if match:
            job_id = match.group(1)
            current_job_id = job_id
            job_running = True
            button.config(state='disabled')
            progressbar.start()
            status_label.config(text=f"{stage} job {job_id} submitted.", fg='blue')
            save_gui_state()
            root.after(200, lambda: check_job_status(job_id, stage))
        else:
            status_label.config(text=f"No job ID found in sbatch output: {stdout}", fg='red')
    except subprocess.CalledProcessError as cpe:
        stderr = cpe.stderr.strip() if cpe.stderr else str(cpe)
        status_label.config(text=f"sbatch error: {stderr}", fg='red')
        button.config(state='normal')  # re-enable on failure
    except Exception as e:
        status_label.config(text=f"Error starting job: {e}", fg='red')
        button.config(state='normal')  # re-enable on failure


# --------------------------
# Check SLURM job status without blocking GUI
# --------------------------
def check_job_status(job_id, stage="QC"):
    global job_running, current_job_id

    button = QC_run if stage=="QC" else filter_run

    try:
        result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
        stdout = result.stdout
    except Exception as e:
        status_label.config(text=f"Error checking job status: {e}", fg='red')
        root.after(5000, lambda: check_job_status(job_id, stage))
        return

    if job_id in stdout:
        status_label.config(text=f"{stage} job {job_id} is still running...", fg='orange')
        root.after(5000, lambda: check_job_status(job_id, stage))
    else:
        status_label.config(text=f"{stage} job {job_id} is complete!", fg='green')
        progressbar.stop()
        progressbar['value'] = 100
        job_running = False
        current_job_id = None
        save_gui_state()
        if stage=="QC":
            next_process()  # enable QC result view buttons
        check_all_ready()   # Re-check if Run QC can be run again
        check_filter_ready() # Re-check if Run Filtering can now be enabled


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
                    f.write("        'results/filtering_done.txt'\n")
                    inside_all = True
                elif inside_all and line.strip().startswith("input:"):
                    continue
                else:
                    f.write(line)
    except Exception as e:
        status_label.config(text=f"Error patching Snakefile: {e}", fg='red')

def enable_filtering_stage():
    for child in scrollable_frame.winfo_children():
        if isinstance(child, ttk.Button) and child.cget("text") == "Save":
            child.config(state='normal')
    patch_snakefile_for_filtering()
    status_label.config(text="Filtering stage enabled.", fg='blue')

# --------------------------
# GUI state persistence
# --------------------------
def save_gui_state():
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
        "current_job_id": current_job_id,
        "job_running": job_running
    }
    try:
        with open(STATE_FILE, "w") as f:
            json.dump(state, f)
    except Exception as e:
        status_label.config(text=f"Error saving GUI state: {e}", fg='red')

def load_gui_state():
    global current_job_id, job_running
    if not os.path.exists(STATE_FILE):
        return
    try:
        with open(STATE_FILE, "r") as f:
            state = json.load(f)
    except Exception as e:
        status_label.config(text=f"Error loading GUI state: {e}", fg='red')
        return

    if 'data_dir_entry' in globals():
        data_dir_entry.delete(0, tk.END)
        data_dir_entry.insert(0, state.get("data_dir", ""))
    if 'metadata_dir_entry' in globals():
        metadata_dir_entry.delete(0, tk.END)
        metadata_dir_entry.insert(0, state.get("metadata_table", ""))
    if 'sample_key_entry' in globals():
        sample_key_entry.delete(0, tk.END)
        sample_key_entry.insert(0, state.get("sample_key", ""))
    if 'seq_batch_entry' in globals():
        seq_batch_entry.delete(0, tk.END)
        seq_batch_entry.insert(0, state.get("seq_batch_key", ""))

    cellranger_saved.set(state.get("cellranger_saved", False))
    metadata_saved.set(state.get("metadata_saved", False))
    sample_key_saved.set(state.get("sample_key_saved", False))
    seq_batch_key_saved.set(state.get("seq_batch_key_saved", False))

    current_job_id = state.get("current_job_id", None)
    job_running = state.get("job_running", False)

    mito_percent_thresh.set(state.get("mito_percent_thresh", ""))
    ribo_percent_thresh.set(state.get("ribo_percent_thresh", ""))
    doublet_thresh.set(state.get("doublet_thresh", ""))
    min_genes_per_cell.set(state.get("min_genes_per_cell", ""))

    check_all_ready()

    if current_job_id and job_running:
        QC_run.config(state='disabled')
        filter_run.config(state='disabled')
        progressbar.start()
        root.after(200, lambda: check_job_status(current_job_id))

root.protocol("WM_DELETE_WINDOW", lambda: (save_gui_state(), root.destroy()))

# --------------------------
# Interface Text / Header
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Single-cell QC Pipeline', style='Header.Label')
header.grid(row=0, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Optional Image (logo)
# --------------------------
try:
    img = Image.open("images/VIRGIL.png")
    img = img.resize((150, 150))
    photo = ImageTk.PhotoImage(img)
    # Keep a reference to PhotoImage to avoid garbage collection (tk quirk)
    logo_label = tk.Label(scrollable_frame, image=photo)
    logo_label.image = photo
    logo_label.grid(row=0, column=4, padx=10, sticky='e')
except Exception:
    # Fail silently if image missing; the GUI still works
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

# Keep references to specific entries by name so load/save can use them
for i, (label_text, var_name, flag_var) in enumerate(entries):
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=i+1, column=0, sticky='e', padx=10, pady=5)
    entry = ttk.Entry(scrollable_frame, width=50)
    entry.grid(row=i+1, column=1, sticky='w')
    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=i+1, column=3, sticky='w')
    # Use default arguments in lambda to capture the current objects
    btn = ttk.Button(scrollable_frame, text='Save', command=lambda v=var_name, e=entry, s=status, f=flag_var: fill(v, e.get(), s, f))
    btn.grid(row=i+1, column=2, sticky='w', padx=5)

    # Store references to the named entries in globals (used by save/load)
    if var_name == "data_dir":
        data_dir_entry = entry
    elif var_name == "metadata_table":
        metadata_dir_entry = entry
    elif var_name == "sample_key":
        sample_key_entry = entry
    elif var_name == "seq_batch_key":
        seq_batch_entry = entry

# --------------------------
# Run + Progress + Status
# --------------------------
QC_run = ttk.Button(scrollable_frame, text='Run QC!', command=start_snakemake_job, state='disabled')
QC_run.grid(row=6, column=3, pady=10, padx=10, sticky='w')

progressbar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
progressbar.grid(row=6, column=0, columnspan=3, padx=10, sticky='ew')

status_label = tk.Label(scrollable_frame, text="Waiting for inputs...", fg="black")
status_label.grid(row=7, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Interface Text / Header
# --------------------------
header = ttk.Label(scrollable_frame, text='scVIRGILS - Filtering', style='Header.Label')
header.grid(row=9, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# QC View buttons (disabled until pipeline finishes)
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
    # Each button opens a local file (png) using the default system viewer
    btn = ttk.Button(
        scrollable_frame,
        text=f"View {text} QC",
        command=lambda p=path, k=qc_key_map[text]: (webbrowser.open(f"file://{os.path.abspath(p)}"), mark_qc_viewed(k)),
        state='disabled'
    )
    btn.grid(row=10+i, column=0, columnspan=2, padx=10, pady=5, sticky='w')

    # Keep names for later enabling
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
# Additional threshold entries (second block)
# --------------------------
base_row = 16  # adjust if needed so it doesn't overlap other widgets

threshold_entries = [
    ("Mitochondria % threshold (e.g. 20)", "mito_percent_thresh", mito_percent_thresh),
    ("Ribosomal % threshold (e.g. 20)", "ribo_percent_thresh", ribo_percent_thresh),
    ("Doublet threshold (e.g. 0.15)", "doublet_thresh", doublet_thresh),
    ("Minimum genes per cell (e.g. 200)", "min_genes_per_cell", min_genes_per_cell)
]

# Store references to Entry widgets for later
entry_widgets = {}

for i, (label_text, var_name, tk_var) in enumerate(threshold_entries):
    row_idx = base_row + i
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=row_idx, column=0, sticky='e', padx=10, pady=5)

    entry = ttk.Entry(scrollable_frame, width=50, textvariable=tk_var)
    entry.grid(row=row_idx, column=1, sticky='w')

    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=row_idx, column=3, sticky='w')

    # Save button disabled initially; will be enabled by Move to Filtering
    btn = ttk.Button(
        scrollable_frame,
        text='Save',
        command=lambda v=var_name, e=entry, s=status: (fill(v, e.get(), s, None, numeric=True), check_filter_ready()),
        state='disabled'
    )
    btn.grid(row=row_idx, column=2, sticky='w', padx=5)

    # Keep references for later use
    entry_widgets[var_name] = entry

# Example: retrieving a threshold later
# Note: calling .get() here right away will return what's currently in the entry (probably empty).
mito_value = entry_widgets["mito_percent_thresh"].get()  # this is fine but will be "" until user types and saves
ribo_value = entry_widgets["ribo_percent_thresh"].get()
doublet_value = entry_widgets["doublet_thresh"].get()
min_genes_per_cell_value = entry_widgets["min_genes_per_cell"].get()

# --------------------------
# Run + Progress + Status
# --------------------------
filter_run = ttk.Button(scrollable_frame, text='Run Filtering!', command=lambda: start_snakemake_job(stage="Filtering"), state='disabled')
filter_run.grid(row=20, column=3, pady=10, padx=10, sticky='w')

progressbar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
progressbar.grid(row=20, column=0, columnspan=3, padx=10, sticky='ew')

status_label = tk.Label(scrollable_frame, text="Waiting for inputs...", fg="black")
status_label.grid(row=20, column=0, columnspan=4, pady=10, sticky='w')

# --------------------------
# Load any previous GUI state after UI creation
# --------------------------
load_gui_state()

# --------------------------
# Start the GUI event loop
# --------------------------
root.mainloop()


