import tkinter as tk
from tkinter import ttk
import subprocess
import re
from PIL import Image, ImageTk
import webbrowser
import json
import os

# Globals
current_job_id = None
job_running = False

# File to store GUI state
STATE_FILE = "gui_state.json"

# Initialize root
root = tk.Tk()
root.geometry('1200x700')
root.minsize(900, 600)
root.title('scVIRGILS')

# Scrollable Canvas Frame
container = tk.Frame(root)
container.pack(fill=tk.BOTH, expand=True)

canvas = tk.Canvas(container, borderwidth=0)
canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
canvas.configure(yscrollcommand=scrollbar.set)

scrollable_frame = tk.Frame(canvas)
scrollable_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")

def _on_mousewheel(event):
    canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

def on_frame_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

def on_canvas_configure(event):
    canvas.itemconfig(scrollable_window, width=event.width)

scrollable_frame.bind("<Configure>", on_frame_configure)
canvas.bind("<Configure>", on_canvas_configure)
root.bind_all("<MouseWheel>", _on_mousewheel)

# Input state vars
cellranger_saved = tk.BooleanVar(value=False)
metadata_saved = tk.BooleanVar(value=False)
sample_key_saved = tk.BooleanVar(value=False)
seq_batch_key_saved = tk.BooleanVar(value=False)

# Style
style = ttk.Style()
style.configure('Header.Label', font=('Helvetica', 16, 'bold'))
style.map('TButton', foreground=[('disabled', '#808080')], background=[('disabled', '#d3d3d3')])

# Fill Snakefile Function
def fill(variable_name, entered_path, status_label, flag_var):
    target_file = "snakefile"
    new_line = f'{variable_name} = "{entered_path}"\n'
    updated = False
    try:
        with open(target_file, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        status_label.config(text=f"Error: File '{target_file}' not found.", foreground='red')
        return

    with open(target_file, 'w') as file:
        for line in lines:
            if line.strip().startswith(f'{variable_name} ='):
                file.write(new_line)
                updated = True
            else:
                file.write(line)
        if not updated:
            file.write(new_line)

    status_label.config(text=f"{variable_name} saved!", fg='green')
    flag_var.set(True)
    check_all_ready()

# Readiness check
def check_all_ready():
    if cellranger_saved.get() and metadata_saved.get() and sample_key_saved.get() and seq_batch_key_saved.get():
        if not job_running:
            QC_run.config(state='normal')
        else:
            QC_run.config(state='disabled')
    else:
        QC_run.config(state='disabled')

# Submit Snakemake job
def start_snakemake_job():
    global current_job_id, job_running
    try:
        result = subprocess.run(['sbatch', 'snakemake.sh'], capture_output=True, check=True, text=True)
        match = re.search(r'(\d+)', result.stdout.strip())
        if match:
            job_id = match.group(1)
            current_job_id = job_id
            job_running = True
            QC_run.config(state='disabled')
            status_label.config(text=f"Submitted job {job_id}", fg='blue')
            progressbar.start()
            save_gui_state()
            check_job_status(job_id)
        else:
            status_label.config(text=f"No job ID found in sbatch output: {result.stdout}", fg='red')
    except Exception as e:
        status_label.config(text=f"Error starting job: {e}", fg='red')

# Check SLURM job status
def check_job_status(job_id):
    global job_running, current_job_id
    result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
    if job_id in result.stdout:
        status_label.config(text=f"Job {job_id} is still running...", fg='orange')
        root.after(5000, lambda: check_job_status(job_id))
    else:
        status_label.config(text=f"Job {job_id} is complete!", fg='green')
        progressbar.stop()
        progressbar['value'] = 100
        job_running = False
        current_job_id = None
        save_gui_state()
        next_process()

# Next step enablement
def next_process():
    open_mito_qc.config(state='normal')
    open_ribo_qc.config(state='normal')
    open_gene_qc.config(state='normal')
    open_doublet_qc.config(state='normal')
    open_genes_by_counts_qc.config(state='normal')

# GUI State persistence
def save_gui_state():
    state = {
        "data_dir": data_dir_entry.get(),
        "metadata_table": metadata_dir_entry.get(),
        "sample_key": sample_key_entry.get(),
        "seq_batch_key": seq_batch_entry.get(),
        "cellranger_saved": cellranger_saved.get(),
        "metadata_saved": metadata_saved.get(),
        "sample_key_saved": sample_key_saved.get(),
        "seq_batch_key_saved": seq_batch_key_saved.get(),
        "mito_perc_threshold": mito_perc_threshold.get(),
        "current_job_id": current_job_id,
        "job_running": job_running
    }
    with open(STATE_FILE, "w") as f:
        json.dump(state, f)

def load_gui_state():
    global current_job_id, job_running
    if not os.path.exists(STATE_FILE):
        return
    with open(STATE_FILE, "r") as f:
        state = json.load(f)
    data_dir_entry.insert(0, state.get("data_dir", ""))
    metadata_dir_entry.insert(0, state.get("metadata_table", ""))
    sample_key_entry.insert(0, state.get("sample_key", ""))
    seq_batch_entry.insert(0, state.get("seq_batch_key", ""))
    cellranger_saved.set(state.get("cellranger_saved", False))
    metadata_saved.set(state.get("metadata_saved", False))
    sample_key_saved.set(state.get("sample_key_saved", False))
    seq_batch_key_saved.set(state.get("seq_batch_key_saved", False))
    current_job_id = state.get("current_job_id", None)
    job_running = state.get("job_running", False)
    mito_perc_threshold = state.get("mito_perc_threshold", False)
    check_all_ready()
    if current_job_id:
        QC_run.config(state='disabled')
        progressbar.start()
        check_job_status(current_job_id)

# Hook close
root.protocol("WM_DELETE_WINDOW", lambda: (save_gui_state(), root.destroy()))

# Interface Text
header = ttk.Label(scrollable_frame, text='scVIRGILS - Single-cell QC Pipeline', style='Header.Label')
header.grid(row=0, column=0, columnspan=4, pady=10, sticky='w')

# Image
try:
    img = Image.open("images/VIRGIL.png")
    img = img.resize((150, 150))
    photo = ImageTk.PhotoImage(img)
    tk.Label(scrollable_frame, image=photo).grid(row=0, column=4, padx=10, sticky='e')
except:
    pass

# Entries and Buttons
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
    if var_name == "data_dir": data_dir_entry = entry
    elif var_name == "metadata_table": metadata_dir_entry = entry
    elif var_name == "sample_key": sample_key_entry = entry
    elif var_name == "seq_batch_key": seq_batch_entry = entry

QC_run = ttk.Button(scrollable_frame, text='Run QC!', command=start_snakemake_job, state='disabled')
QC_run.grid(row=6, column=3, pady=10, padx=10, sticky='w')

progressbar = ttk.Progressbar(scrollable_frame, mode='indeterminate')
progressbar.grid(row=6, column=0, columnspan=3, padx=10, sticky='ew')

status_label = tk.Label(scrollable_frame, text="Waiting for inputs...", fg="black")
status_label.grid(row=7, column=0, columnspan=4, pady=10, sticky='w')

# QC View buttons
view_buttons = [
    ("% Mitochondria", 'figures/QC_mito_pct.png'),
    ("% Ribosomal", 'figures/QC_ribo_pct.png'),
    ("Gene Counts", 'figures/QC_gene_counts.png'),
    ("Doublet Score", 'figures/QC_doublet.png'),
    ("Genes by Counts", 'figures/QC_genes_by_counts.png')
]

for i, (text, path) in enumerate(view_buttons):
    btn = ttk.Button(scrollable_frame, text=f"View {text} QC", command=lambda p=path: webbrowser.open(f"file://{os.path.abspath(p)}"), state='disabled')
    btn.grid(row=9+i, column=0, columnspan=2, padx=10, pady=5, sticky='w')
    if i == 0: open_mito_qc = btn
    elif i == 1: open_ribo_qc = btn
    elif i == 2: open_gene_qc = btn
    elif i == 3: open_doublet_qc = btn
    elif i == 4: open_genes_by_counts_qc = btn

# Add Filtering thresholds for filtering
# Entries and Buttons
# Define the entries you want as tuples (label, var_name)
next_entries = [
    ("Mitochondria % threshold", "mito_percent_thresh"),
    ("Ribosomal % threshold", "ribo_percent_thresh"),
    # add more entries here
]

# Dictionary to hold references to entry widgets
entry_widgets = {}

for i, (label_text, var_name) in enumerate(next_entries):
    ttk.Label(scrollable_frame, text=label_text, wraplength=400).grid(row=i+1, column=0, sticky='e', padx=10, pady=5)
    
    entry = ttk.Entry(scrollable_frame, width=50)
    entry.grid(row=i+1, column=1, sticky='w')
    
    status = tk.Label(scrollable_frame, text="", anchor='w')
    status.grid(row=i+1, column=3, sticky='w')
    
    btn = ttk.Button(scrollable_frame, text='Save',
                     command=lambda v=var_name, e=entry, s=status: fill(v, e.get(), s))
    btn.grid(row=i+1, column=2, sticky='w', padx=5)
    
    # Save entry widget for later use if needed
    entry_widgets[var_name] = entry

# Now, for example, to get the Mito threshold value somewhere else, use:
mito_value = entry_widgets["mito_percent_thresh"].get()




# Load previous state
load_gui_state()

# Mainloop
root.mainloop()

