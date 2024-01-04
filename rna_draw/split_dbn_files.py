import os
import shutil
import click
import numpy as np

def common_replaces(job_text, args):
    for k, v in args.items():
        job_text = job_text.replace('{' + k + '}', str(v))
    return job_text

@click.group()
def cli():
    pass

@cli.command()  # @cli, not @click!
@click.argument('cmd_file', type=click.Path(exists=True))
@click.argument('dbn_dir', type=click.Path(exists=True))
@click.option('--n', 'n_splits', default=1000, help='Number of splits')
@click.option('-rd', '--run-dir', default='sortedDBNFiles', help='where to store the sorted DBN files')
@click.option('--time', 'JOB_TIME', default="24:00:00", help='how long should the jobs be')
def split_dbn(cmd_file, dbn_dir, n_splits, run_dir, **args):
    run_dir = '/work/yesselmanlab/nklein' #os.path.abspath(run_dir)
    job_text = open(cmd_file).read()
    args['WORK_DIR'] = os.path.abspath('.')
    job_text = common_replaces(job_text, args)
    dbn_files = os.listdir(dbn_dir)
    dbn_files = [f for f in dbn_files if f.endswith('.dbn')]
    dbn_file_splits = np.array_split(dbn_files, n_splits)
    f_sum = open('README_SUBMIT', 'w')
    os.makedirs(run_dir, exist_ok=True)

    for i, file_group in enumerate(dbn_file_splits):
        job_dir = os.path.join(run_dir, str(i))

        success_dir = os.path.join(job_dir, "Success")
        fails_dir = os.path.join(job_dir, "Fails")
        os.makedirs(success_dir, exist_ok=True)
        os.makedirs(fails_dir, exist_ok=True)

        os.makedirs(job_dir, exist_ok=True)
        for file in file_group:
            shutil.copy(os.path.join(dbn_dir, file), job_dir)
        submit_path = os.path.join(job_dir, 'job.sh')
        f = open(submit_path, 'w')
        job_text_final = job_text.replace('{JOB_DIR}', job_dir)
        f.write(job_text_final)
        f.close()
        f_sum.write(f"sbatch {submit_path}\n")
    f_sum.close()

if __name__ == '__main__':
    cli()
