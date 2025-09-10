#!/usr/bin/env python3
"""
scVIRGILS - Tkinter GUI for running a Snakemake QC pipeline
This file is a debugged, heavily commented version of the user's original script.
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
current_job_id = None      # ID of a currently-submitted SLURM job (string)
job_running = False       # Bool - is a job known to be running?
STATE_FILE = "gui_state.json"

# --------------------------
# Basic Tk root + layout
# --------------------------
root = tk.Tk()
root.geometry('1200x700')
root.minsize(900, 600)
root.title('scVIRGILS')

# Create a scrollable canvas frame pattern:
# container -> canvas -> scrollable_frame (a frame inside canvas)
container = tk.Frame(root)
container.pack(fill=tk.BOTH, expand=True)

canvas = tk.Canvas(container, borderwidth=0)
canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
canvas.configure(yscrollcommand=scrollbar.set)

# The actual frame where widgets live
scrollable_frame = tk.Frame(canvas)
scrollable_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

# ------------ scrolling handlers ------------
def _on_mousewheel(event):
    """
    Generic mousewheel scrolling. Note: event.delta behaves differently on Mac/Linux.
    This handler keeps the behaviour similar to Windows; you can add platform checks
    if you want separate handling for Linux (Button-4/5) or Mac.
    """
    canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

def on_frame_configure(event):
    """Update the scrollregion when the scrollable_frame changes size."""
    canvas.configure(scrollregion=canvas.bbox("all"))

def on_canvas_configure(event):
    """Make the inner window the same width as the canvas (responsive)."""
    canvas.itemconfig(scrollable_window, width=event.width)

scrollable_frame.bind("<Configure>", on_frame_configure)
canvas.bind("<Configure>", on_canvas_configure)
root.bind_all("<MouseWheel>", _on_mousewheel)

# --------------------------
# Input state variables
# --------------------------
# These booleans track whether the user has saved the 4 required inputs.
cellranger_saved = tk.BooleanVar(value=False)
metadata_saved = tk.BooleanVar(value=False)
sample_key_saved = tk.BooleanVar(value=False)
seq_batch_key_saved = tk.BooleanVar(value=False)

# An example threshold variable referenced by save/load functions (must exist)
mito_perc_threshold = tk.StringVar(value="")  # was missing in original code

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
def fill(variable_name, entered_path, status_label=None, flag_var=None):
    """
    Write or update a line of the form:
        variable_name = "entered_path"
    in a file called 'snakefile' (same directory). If the variable is not present,
    append it.

    Parameters:
    - variable_name (str): variable to set inside the snakefile
    - entered_path (str): the string value to write
    - status_label (tk.Label or None): UI label to update with messages
    - flag_var (tk.BooleanVar or None): optional boolean Tk variable to set True on success
    """
    target_file = "snakefile"  # change to "Snakefile" if that's the real filename
    # Represent the assignment exactly as in Snakefile e.g. foo = "bar"
    new_line = f'{variable_name} = "{entered_path}"\n'
    updated = False

    try:
        with open(target_file, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        # If the Snakefile doesn't exist, inform the user via status_label (if provided)
        if status_label is not None:
            status_label.config(text=f"Error: File '{target_file}' not found.", fg='red')
        return

    # Write back updating or appending the variable assignment
    try:
        with open(target_file, 'w') as file:
            for line in lines:
                # Match the variable assignment at start of a stripped line
                if line.strip().startswith(f'{variable_name} ='):
                    file.write(new_line)
                    updated = True
                else:
                    file.write(line)
            if not updated:
                # If variable wasn't in file, append it at the end
                file.write(new_line)
    except Exception as e:
        if status_label is not None:
            status_label.config(text=f"Error writing '{target_file}': {e}", fg='red')
        return

    # Update UI + flag if present
    if status_label is not None:
        status_label.config(text=f"{variable_name} saved!", fg='green')
    if flag_var is not None:
        try:
            flag_var.set(True)
        except Exception:
            # If flag_var isn't a Tk variable, ignore silently (or log)
            pass

    # Reevaluate whether all required inputs are present so we can enable "Run QC!"
    check_all_ready()

# --------------------------
# Enable/Disable Run button based on readiness
# --------------------------
def check_all_ready():
    """
    Enable QC_run only when all four input flags are True and no job is running.
    Note: job_running is a global boolean that is updated by job submission / status checks.
    """
    global job_running
    if cellranger_saved.get() and metadata_saved.get() and sample_key_saved.get() and seq_batch_key_saved.get():
        if not job_running:
            QC_run.config(state='normal')
        else:
            QC_run.config(state='disabled')
    else:
        QC_run.config(state='disabled')

# --------------------------
# Start a Snakemake job via sbatch
# --------------------------
def start_snakemake_job():
    """
    Submits the 'snakemake.sh' script using sbatch. Parses the returned job ID and
    starts periodic status checks. All subprocess calls deliberately run with text=True
    so output is captured as strings.

    Note: sbatch stdout typically contains a job id like "Submitted batch job 12345".
    We attempt to extract digits from stdout.
    """
    global current_job_id, job_running

    try:
        result = subprocess.run(['sbatch', 'snakemake.sh'], capture_output=True, check=True, text=True)
        stdout = result.stdout.strip()
        # Extract the first contiguous string of digits as the job id
        match = re.search(r'(\d+)', stdout)
        if match:
            job_id = match.group(1)
            current_job_id = job_id
            job_running = True
            QC_run.config(state='disabled')
            status_label.config(text=f"Submitted job {job_id}", fg='blue')
            progressbar.start()
            save_gui_state()
            # Schedule the first check shortly after returning to the event loop to avoid blocking
            root.after(200, lambda: check_job_status(job_id))
        else:
            status_label.config(text=f"No job ID found in sbatch output: {stdout}", fg='red')
    except subprocess.CalledProcessError as cpe:
        # If sbatch returns a non-zero exit code, show stderr for diagnosis
        stderr = cpe.stderr.strip() if cpe.stderr else str(cpe)
        status_label.config(text=f"sbatch error: {stderr}", fg='red')
    except Exception as e:
        status_label.config(text=f"Error starting job: {e}", fg='red')

# --------------------------
# Check SLURM job status without blocking GUI
# --------------------------
def check_job_status(job_id):
    """
    Uses squeue -j <job_id> to check for job presence. If squeue returns the job id,
    reschedule another check. If the job is not in squeue output, we assume it completed.

    Note: This function performs subprocess.run(...) which is a blocking call but is
    executed from the Tk event loop via root.after to avoid freezing the UI for long.
    """
    global job_running, current_job_id

    try:
        # This call is short but could block if Slurm unresponsive; keeping it simple.
        result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
    except Exception as e:
        status_label.config(text=f"Error checking job status: {e}", fg='red')
        # try again later
        root.after(5000, lambda: check_job_status(job_id))
        return

    stdout = result.stdout
    if job_id in stdout:
        # Still running -> check again in 5 seconds
        status_label.config(text=f"Job {job_id} is still running...", fg='orange')
        root.after(5000, lambda: check_job_status(job_id))
    else:
        # Job no longer in queue; mark complete
        status_label.config(text=f"Job {job_id} is complete!", fg='green')
        progressbar.stop()
        progressbar['value'] = 100
        job_running = False
        current_job_id = None
        save_gui_state()
        next_process()

# --------------------------
# Enable the next set of QC buttons after run completes
# --------------------------
def next_process():
    """Enable the QC result view buttons when pipeline is finished."""
    open_mito_qc.config(state='normal')
    open_ribo_qc.config(state='normal')
    open_gene_qc.config(state='normal')
    open_doublet_qc.config(state='normal')
    open_genes_by_counts_qc.config(state='normal')

# --------------------------
# GUI State persistence
# --------------------------
def save_gui_state():
    """
    Save selected GUI values into a JSON file so GUI state can be restored later.
    Save only serializable values (strings, booleans). Do not attempt to save Tk objects.
    """
    state = {
        "data_dir": data_dir_entry.get() if 'data_dir_entry' in globals() else "",
        "metadata_table": metadata_dir_entry.get() if 'metadata_dir_entry' in globals() else "",
        "sample_key": sample_key_entry.get() if 'sample_key_entry' in globals() else "",
        "seq_batch_key": seq_batch_entry.get() if 'seq_batch_entry' in globals() else "",
        "cellranger_saved": cellranger_saved.get(),
        "metadata_saved": metadata_saved.get(),
        "sample_key_saved": sample_key_saved.get(),
        "seq_batch_key_saved": seq_batch_key_saved.get(),
        "mito_perc_threshold": mito_perc_threshold.get(),
        "current_job_id": current_job_id,
        "job_running": job_running
    }
    try:
        with open(STATE_FILE, "w") as f:
            json.dump(state, f)
    except Exception as e:
        # If saving fails, surface to status_label but continue
        status_label.config(text=f"Error saving GUI state: {e}", fg='red')

def load_gui_state():
    """
    Load GUI state from disk (if present). This should be called after all Entry widgets
    and corresponding Tk variables have been created, otherwise insert/set calls will fail.
    """
    global current_job_id, job_running

    if not os.path.exists(STATE_FILE):
        return

    try:
        with open(STATE_FILE, "r") as f:
            state = json.load(f)
    except Exception as e:
        status_label.config(text=f"Error loading GUI state: {e}", fg='red')
        return

    # Populate entry widgets only if they exist
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

    # Restore boolean flags
    cellranger_saved.set(state.get("cellranger_saved", False))
    metadata_saved.set(state.get("metadata_saved", False))
    sample_key_saved.set(state.get("sample_key_saved", False))
    seq_batch_key_saved.set(state.get("seq_batch_key_saved", False))

    # Restore job metadata
    current_job_id = state.get("current_job_id", None)
    job_running = state.get("job_running", False)

    # Set the Tk var for mito threshold rather than overwriting the object
    mito_perc_threshold.set(state.get("mito_perc_threshold", ""))

    # Re-check enabling logic for the Run button
    check_all_ready()

    # If a job was running when we saved, resume checking status (non-blocking)
    if current_job_id and job_running:
        QC_run.config(state='disabled')
        progressbar.start()
        root.after(200, lambda: check_job_status(current_job_id))

# When the window is closed, save state then destroy window
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
# QC View buttons (disabled until pipeline finishes)
# --------------------------
view_buttons = [
    ("% Mitochondria", 'figures/QC_mito_pct.png'),
    ("% Ribosomal", 'figures/QC_ribo_pct.png'),
    ("Gene Counts", 'figures/QC_gene_counts.png'),
    ("Doublet Score", 'figures/QC_doublet.png'),
    ("Genes by Counts", 'figures/QC_genes_by_counts.png')
]

for i, (text, path) in enumerate(view_buttons):
    # Each button opens a local file (png) using the default system viewer
    btn = ttk.Button(scrollable_frame,
                     text=f"View {text} QC",
                     command=lambda p=path: webbrowser.open(f"file://{os.path.abspath(p)}"),
                     state='disabled')
    btn.grid(row=9+i, column=0, columnspan=2, padx=10, pady=5, sticky='w')

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

# --------------------------
# Additional threshold entries (second block)
# --------------------------
# Start these rows well below the first block to avoid overlapping widget grid indices.
base_row = 15

next_entries = [
    ("Mitochondria % threshold", "mito_percent_thresh"),
    ("Ribosomal % threshold", "ribo_percent_thresh"),
    # add more entries here
]

# Keep references to these entry widgets if you need their values later
entry_widgets = {}

for i, (label_text, var_name) in enumerate(next_entries):
    row_idx = base_row + i
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=row_idx, column=0, sticky='e', padx=10, pady=5)

    entry = ttk.Entry(scrollable_frame, width=50)
    entry.grid(row=row_idx, column=1, sticky='w')

    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=row_idx, column=3, sticky='w')

    # Use a version of fill that does not require a flag_var (it is optional now).
    btn = ttk.Button(scrollable_frame, text='Save', command=lambda v=var_name, e=entry, s=status: fill(v, e.get(), s, None))
    btn.grid(row=row_idx, column=2, sticky='w', padx=5)

    # Save reference for later usage
    entry_widgets[var_name] = entry

# Example: retrieving a threshold later
# Note: calling .get() here right away will return what's currently in the entry (probably empty).
mito_value = entry_widgets["mito_percent_thresh"].get()  # this is fine but will be "" until user types and saves

# --------------------------
# Load any previous GUI state after UI creation
# --------------------------
load_gui_state()

# --------------------------
# Start the GUI event loop
# --------------------------
root.mainloop()


