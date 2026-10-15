import os
import signal
import subprocess


def kill_job(job):
    if not hasattr(job, 'pid'):
        return
    try:
        os.killpg(os.getpgid(job.pid), signal.SIGTERM)
    except ProcessLookupError:
        pass


def execute_script_job(log_path, script_path, active_jobs, name):
    with open(log_path, 'w') as lf:
        p = subprocess.Popen(["bash", script_path], stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
        active_jobs[name] = p
        p.wait()
        active_jobs.pop(name, None)
