import glob
import os
import subprocess
import typing

from . import logger as log
from .configuration import config
from .functions import substitute_template


def create_bash_files(files: list, command: str) -> list:
    """
    Generate submit scripts for each input file files and for calculation on a
    cluster.

    Parameters
    ----------
    files : list of str
        List of input files on which the command will be executed.
    command : str
        Command to be executed and used for submitting scripts to a cluster.

    Returns
    -------
    list of str
        List containing paths to the generated bash scripts.
    """
    bash_files = []
    basepath = os.getcwd()
    command = substitute_template(command, config)[0]

    for f in files:
        fpath = os.path.dirname(f) or "."
        os.chdir(fpath)
        try:
            infilename = os.path.basename(f)
            cmd = command.strip() + " " + infilename
            log.info(f"Now building script for {infilename}: {f}", "Bash Wrapper")
            result = subprocess.run(cmd, shell=True, check=False)
            if getattr(result, "returncode", 0) != 0:
                log.error(
                    f"Command failed for {infilename} (exit code "
                    f"{result.returncode}): {cmd}",
                    "Bash Wrapper",
                )
            log.debug(f"Executing command : {cmd}", "Bash Wrapper")

            fname_no_ext = os.path.splitext(infilename)[0]
            all_shs = glob.glob("*.sh")
            all_shs.extend(glob.glob("*.sbatch"))

            local_execs = [os.path.abspath(sh) for sh in all_shs if fname_no_ext in sh]
            bash_files.extend(local_execs)
        finally:
            os.chdir(basepath)
    return bash_files


def make_send_script(bashfiles: list, send_command: str, sendscript: typing.IO):
    """
    Generate a script for sending all jobscripts to a cluster.

    Parameters
    ----------
    bashfiles : list of str
        List of paths to the bash scripts to be sent and executed.
    send_command : str
        Command used for sending and executing bash scripts on a cluster.
    sendscript : typing.IO
        File-like object representing the script file to be generated.
    """
    send_command = substitute_template(send_command, config)[0]
    sendscript_content = (
        "#!/bin/bash\n"
        "function cd_and_sbatch() {\n"
        '    local script_file="$1"\n'
        '    local folder="$2"\n'
        '    echo "Sending $script_file"\n'
        '    cd "$folder"\n'
        f'    {send_command.strip()} "$script_file"\n'
        "}\n\n"
    )

    for f in bashfiles:
        fpath = os.path.abspath(os.path.dirname(f))
        infilename = os.path.basename(f)

        addendum = f'cd_and_sbatch "{infilename}" "{fpath}"\n'
        sendscript_content += addendum

    sendscript.write(sendscript_content)
