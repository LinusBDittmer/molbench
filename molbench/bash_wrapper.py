import os
import glob
import subprocess
import typing
from . import logger as log
from .functions import substitute_template
from . import config


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

    Notes
    -----
    This function generates bash scripts that submit them to a cluster for
    execution.

    Explanation
    -----------
    - `bash_files`: List to store paths of generated bash scripts.
    - `basepath`: Get the current working directory.
    - `command`: Command used for submitting scripts to a cluster.
    - Iterate over each input file:
        - `fpath`: Get the directory path of the current file.
        - Change the current working directory to `fpath`.
        - `infilename`: Extract the base name of the current file.
        - Construct the command to be executed on the current file.
        - Log a message indicating the script generation process.
        - Execute the command using subprocess to submit the script to the
          cluster.
        - Log debug information about the executed command.
        - Extract the filename without extension.
        - Find all '.sh' and '.sbatch' files in the current directory.
        - Filter out files related to the current input file.
        - Append absolute paths of matching files to `local_execs`.
        - Extend `bash_files` with the paths of matching files.
        - Change the current working directory back to the base path.
    - Return the list of generated bash script paths.

    """
    bash_files = []
    basepath = os.getcwd()
    command = substitute_template(command, config)[0]

    for f in files:
        fpath = os.path.dirname(f)
        os.chdir(fpath)
        infilename = os.path.basename(f)
        cmd = command.strip() + " " + infilename
        log.info(f"Now building script for {infilename}: {f}", "Bash Wrapper")
        subprocess.run(cmd, shell=True)
        log.debug(f"Executing command : {cmd}", "Bash Wrapper")

        fname_no_ext = os.path.splitext(infilename)[0]
        all_shs = glob.glob("*.sh")
        all_shs.extend(glob.glob("*.sbatch"))

        local_execs = [os.path.abspath(sh) for sh in all_shs
                       if fname_no_ext in sh]
        bash_files.extend(local_execs)
        os.chdir(basepath)
    return bash_files


def make_send_script(bashfiles: list, send_command: str,
                     sendscript: typing.IO):
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

    Notes
    -----
    This function generates a script for sending all jobscripts to a cluster.

    Explanation
    -----------
    - `send_command`: Command used for sending and executing bash scripts on a
      cluster.
    - `sendscript_content`: Initialize the content of the send script with
                            shebang and a function definition.
    - Iterate over each bash script file:
        - `fpath`: Get the absolute directory path of the current bash script.
        - `infilename`: Extract the base name of the current bash script.
        - Construct the addendum to the send script content for sending and
          executing the current bash script.
        - Append the addendum to the sendscript content.
    - Write the sendscript content to the sendscript file.

    """
    send_command = substitute_template(send_command, config)[0]
    sendscript_content = (
        "#!/bin/bash\n"
        "function cd_and_sbatch() {\n"
        "    local script_file=\"$1\"\n"
        "    local folder=\"$2\"\n"
        "    echo \"Sending $script_file\"\n"
        "    cd \"$folder\"\n"
        f"    {send_command.strip()} \"$script_file\"\n"
        "}\n\n"
    )

    for f in bashfiles:
        fpath = os.path.abspath(os.path.dirname(f))
        infilename = os.path.basename(f)

        addendum = f"cd_and_sbatch \"{infilename}\" \"{fpath}\"\n"
        sendscript_content += addendum

    sendscript.write(sendscript_content)
