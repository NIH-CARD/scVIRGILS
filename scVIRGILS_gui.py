import tkinter as tk
from tkinter import ttk
import subprocess
import re
from PIL import Image, ImageTk

# This function helps fill in entered data to the snakefile
def fill(variable_name, entered_path, status_label, flag_var):
    target_file = "scVIRGILS/snakefile"
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

    status_label.config(text=f"Saved: {new_line.strip()}", fg='green')
    flag_var.set(True)
    check_all_ready()

# This is a checker function to make sure the buttons have been clicked and Entries have been saved
def check_all_ready():
    if cellranger_saved.get() and metadata_saved.get() and sample_key_saved.get() and seq_batch_key_saved.get():
        QC_run.config(state='normal')

# This function starts the snakemake job and monitors it
def start_snakemake_job():
    try:
        result = subprocess.run(['sbatch', 'snakemake.sh'], capture_output=True, check=True)
        match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if match:
            job_id = match.group(1)
            status_label.config(text=f"Submitted job {job_id}", fg='blue')
            progressbar.start()
            check_job_status(job_id)
        else:
            status_label.config(text="Failed to parse job ID", fg='red')
    except subprocess.CalledProcessError as e:
        status_label.config(text=f"sbatch failed: {e.stderr}", fg='red')

# This function checks the SLURM job status
def check_job_status(job_id):
    result = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True)
    if job_id in result.stdout:
        status_label.config(text=f"Job {job_id} is still running:", fg='orange')
        root.after(5000, lambda: check_job_status(job_id)) # checks in 5 sec
    else:
        status_label.config(text=f"Job {job_id} is complete!", fg='green')
        progressbar.stop()
        progressbar['value']=100

# ==== GUI ====
root = tk.Tk()
root.geometry('1000x500')
root.title('scVIRGILS')

# Add image to top
for i in range(4):
    root.columnconfigure(i, weight=1)

image_path = "scVIRGILS/images/VIRGIL.png"
image = Image.open(image_path)
image = image.resize((100, 100))
photo = ImageTk.PhotoImage(image)

image_label = tk.Label(root, image=photo)
image_label.image = photo
image_label.grid(row=0, column=4, columnspan=4, pady=(10, 0), sticky='n')

# Saving necessary variables for the succesful entry of information
cellranger_saved = tk.BooleanVar(value=False)
metadata_saved = tk.BooleanVar(value=False)
sample_key_saved = tk.BooleanVar(value=False)
seq_batch_key_saved = tk.BooleanVar(value=False)

# Correct label widget (capital L)
interface_summary = ttk.Label(
    root,
    text='This GUI interface is designed for the execution of scVIRGILS. Each step must be completed before continuing to the next. Good Luck!',
    justify='left'
)
interface_summary.grid(row=0, column=0, columnspan=4, pady=10, padx=10, sticky='w')

# Label for path (CELLRANGER)
path_label_cellranger = ttk.Label(root, text='Path to CELLRANGER files (raw_feature_bc_matrix.h5):')
path_label_cellranger.grid(row=1, column=0, padx=10, pady=10, sticky='e')

# Entry for path (CELLRANGER)
data_dir_entry = ttk.Entry(root)
data_dir_entry.grid(row=1, column=1, padx=10, pady=10, sticky='w')
data_dir_entry.focus()  # This sets the focus of the entry bar

# Status label for feedback (CELLRANGER)
status_label_cellranger = tk.Label(root, text="", anchor='w', justify='left')
status_label_cellranger.grid(row=1, column=3, padx=10, pady=10, sticky='w')

# Path button (CELLRANGER)
path_button_cellranger = ttk.Button(root, text='Save Path', command=lambda: fill("data_dir", data_dir_entry.get(), status_label_cellranger, cellranger_saved))
path_button_cellranger.grid(row=1, column=2, padx=10, pady=10, sticky='w')

# Label for path (METADATA)
path_label_metadata = ttk.Label(root, text='Path the METADATA files (.csv):')
path_label_metadata.grid(row=2, column=0, padx=10, pady=10, sticky='e')

# Entry for path (METADATA)
metadata_dir_entry = ttk.Entry(root)
metadata_dir_entry.grid(row=2, column=1, padx=10, pady=10, sticky='w')

# Status label (METADATA)
status_label_metadata = tk.Label(root, text="", anchor='w', justify='left')
status_label_metadata.grid(row=2, column=3, padx=10, pady=10, sticky='w')

# Path button (METADATA)
path_button_metadata = ttk.Button(root, text='Save Path', command=lambda: fill("metadata_table", metadata_dir_entry.get(), status_label_metadata, metadata_saved))
path_button_metadata.grid(row=2, column=2, padx=10, pady=10, sticky='w')

# Label for SAMPLE_KEY
sample_key_label = ttk.Label(root, text='Key for sample names, required for aggregating while preserving sample info. This must be entered as it appears in metadata (e.g. sample_id):', wraplength=400)
sample_key_label.grid(row=3, column=0, padx=10, ipady=1, sticky='e')

# Entry for SAMPLE_KEY
sample_key_entry = ttk.Entry(root)
sample_key_entry.grid(row=3, column=1, padx=10, pady=10, sticky='w')

# Status label (SAMPLE_KEY)
status_label_sample_key = tk.Label(root, text="", anchor='w', justify='left')
status_label_sample_key.grid(row=3, column=3, padx=10, pady=10, sticky='w')

# Button or SAMPLE_KEY
sample_key_button = ttk.Button(root, text='Save sample_key', command=lambda: fill("sample_key", sample_key_entry.get(), status_label_sample_key, sample_key_saved))
sample_key_button.grid(row=3, column=2, padx=10, pady=10, sticky='w')

# Label for SEQ_BATCH_KEY
seq_batch_key_label = ttk.Label(root, text='Key for sequencing batch, requried for batch correction. This must be entered as it appears in the metadata (e.g. sequencing_round):', wraplength=400)
seq_batch_key_label.grid(row=4, column=0, padx=10, ipady=1, sticky='e')

# Entry for SEQ_BATCH_KEY
seq_batch_entry = ttk.Entry(root)
seq_batch_entry.grid(row=4, column=1, padx=10, pady=10, sticky='w')

# Status label (SEQ_BATCH_KEY)
status_label_seq_batch_key = tk.Label(root, text="", anchor='w', justify='left')
status_label_seq_batch_key.grid(row=4, column=3, padx=10, pady=10, sticky='w')

# Button for SEQ_BATCH_KEY
seq_batch_key_button = ttk.Button(root, text='Save seq_batch_key', command=lambda: fill("seq_batch_key", seq_batch_entry.get(), status_label_seq_batch_key, seq_batch_key_saved))
seq_batch_key_button.grid(row=4, column=2, padx=10, pady=10, sticky='w')

# Button to Start running the snakefile
QC_run = ttk.Button(root, text='STEP 1: Run QC!', command=start_snakemake_job, state='disabled')
QC_run.grid(row=5, column=3, padx=10, pady=10, sticky='w')

# Progress Bar and status label
progressbar = ttk.Progressbar(root, mode='indeterminate', style='Striped.Horizontal.TProgressbar')
progressbar.grid(row=6, column=0, columnspan=4, padx=10, sticky='w')

status_label = tk.Label(root, text="", fg="black", anchor='w')
status_label.grid(row=7, column=0, pady=10, padx=10, sticky='w')


root.mainloop()

