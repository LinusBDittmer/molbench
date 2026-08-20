import io
import os
import subprocess
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_fake_run(sh_suffix=".sh"):
    """Returns a subprocess.run mock that creates a script file in the
    current working directory (which is the input file's directory during
    the create_bash_files call), simulating a submit-script generator."""

    def fake_run(cmd, shell=False, check=False):
        infilename = cmd.strip().split()[-1]
        stem = Path(infilename).stem
        Path(stem + sh_suffix).write_text("#!/bin/bash\n#SBATCH " + stem)

    return fake_run


# ---------------------------------------------------------------------------
# create_bash_files
# ---------------------------------------------------------------------------


def test_create_bash_files_basic(tmp_path, monkeypatch):
    infile = tmp_path / "molA_HF_cc-pvdz.in"
    infile.write_text("fake input")
    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", _make_fake_run(".sh"))
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([str(infile)], "fakegen")
    assert len(result) == 1
    assert result[0].endswith(".sh")
    assert "molA_HF_cc-pvdz" in result[0]


def test_create_bash_files_sbatch_suffix(tmp_path, monkeypatch):
    infile = tmp_path / "mol.in"
    infile.write_text("")
    monkeypatch.setattr(
        "molbench.bash_wrapper.subprocess.run", _make_fake_run(".sbatch")
    )
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([str(infile)], "fakegen")
    assert len(result) == 1
    assert result[0].endswith(".sbatch")


def test_create_bash_files_both_sh_and_sbatch(tmp_path, monkeypatch):
    infile = tmp_path / "mol.in"
    infile.write_text("")

    def fake_run_both(cmd, shell=False, check=False):
        stem = Path(cmd.strip().split()[-1]).stem
        Path(stem + ".sh").write_text("#!/bin/bash")
        Path(stem + ".sbatch").write_text("#!/bin/bash")

    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", fake_run_both)
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([str(infile)], "fakegen")
    assert len(result) == 2
    extensions = {Path(r).suffix for r in result}
    assert ".sh" in extensions
    assert ".sbatch" in extensions


def test_create_bash_files_stem_filter(tmp_path, monkeypatch):
    """A .sh file for a different stem must not be returned."""
    infile = tmp_path / "target.in"
    infile.write_text("")

    def fake_run_creates_both(cmd, shell=False, check=False):
        # Creates target.sh AND unrelated_other.sh
        Path("target.sh").write_text("#!/bin/bash")
        Path("unrelated_other.sh").write_text("#!/bin/bash")

    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", fake_run_creates_both)
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([str(infile)], "fakegen")
    assert all("target" in r for r in result)
    assert not any("unrelated_other" in r for r in result)


def test_create_bash_files_multiple_inputs(tmp_path, monkeypatch):
    files = []
    for name in ("a", "b", "c"):
        f = tmp_path / f"{name}.in"
        f.write_text("")
        files.append(str(f))
    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", _make_fake_run(".sh"))
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files(files, "fakegen")
    assert len(result) == 3


def test_create_bash_files_empty_list(monkeypatch):
    called = []
    monkeypatch.setattr(
        "molbench.bash_wrapper.subprocess.run", lambda *a, **kw: called.append(1)
    )
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([], "fakegen")
    assert result == []
    assert called == []


def test_create_bash_files_restores_cwd(tmp_path, monkeypatch):
    infile = tmp_path / "mol.in"
    infile.write_text("")
    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", lambda *a, **kw: None)
    original_cwd = os.getcwd()
    from molbench.bash_wrapper import create_bash_files

    create_bash_files([str(infile)], "fakegen")
    assert os.getcwd() == original_cwd


def test_create_bash_files_returns_absolute_paths(tmp_path, monkeypatch):
    infile = tmp_path / "mol.in"
    infile.write_text("")
    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", _make_fake_run(".sh"))
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files([str(infile)], "fakegen")
    for path in result:
        assert os.path.isabs(path)


def test_create_bash_files_config_substitution(tmp_path, monkeypatch):
    """[[threads]] in command is resolved from config before subprocess call."""
    infile = tmp_path / "mol.in"
    infile.write_text("")
    received_cmds = []

    def capturing_run(cmd, shell=False, check=False):
        received_cmds.append(cmd)

    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", capturing_run)
    from molbench.bash_wrapper import create_bash_files

    create_bash_files([str(infile)], "run -n [[threads]] mol")
    assert len(received_cmds) == 1
    # [[threads]] should be replaced; the default config has threads=1
    assert "[[threads]]" not in received_cmds[0]
    assert "1" in received_cmds[0]


# ---------------------------------------------------------------------------
# make_send_script
# ---------------------------------------------------------------------------


def test_make_send_script_shebang():
    from molbench.bash_wrapper import make_send_script

    buf = io.StringIO()
    make_send_script([], "sbatch", buf)
    assert buf.getvalue().startswith("#!/bin/bash")


def test_make_send_script_function_definition():
    from molbench.bash_wrapper import make_send_script

    buf = io.StringIO()
    make_send_script([], "sbatch", buf)
    assert "function cd_and_sbatch()" in buf.getvalue()


def test_make_send_script_contains_all_files(tmp_path):
    from molbench.bash_wrapper import make_send_script

    files = [str(tmp_path / "a.sh"), str(tmp_path / "b.sbatch")]
    buf = io.StringIO()
    make_send_script(files, "sbatch", buf)
    content = buf.getvalue()
    assert "a.sh" in content
    assert "b.sbatch" in content


def test_make_send_script_empty_list():
    from molbench.bash_wrapper import make_send_script

    buf = io.StringIO()
    make_send_script([], "sbatch", buf)
    content = buf.getvalue()
    # Only header + function; no cd_and_sbatch calls
    call_lines = [l for l in content.splitlines() if l.startswith("cd_and_sbatch")]
    assert call_lines == []


def test_make_send_script_writes_to_stringio(tmp_path):
    from molbench.bash_wrapper import make_send_script

    buf = io.StringIO()
    make_send_script([str(tmp_path / "x.sh")], "sbatch", buf)
    assert len(buf.getvalue()) > 0


def test_make_send_script_writes_to_real_file(tmp_path):
    from molbench.bash_wrapper import make_send_script

    out = tmp_path / "send.sh"
    with open(out, "w") as f:
        make_send_script([str(tmp_path / "x.sh")], "sbatch", f)
    assert out.stat().st_size > 0


def test_make_send_script_send_command_substitution():
    """[[threads]] in send_command is resolved from config."""
    from molbench.bash_wrapper import make_send_script

    buf = io.StringIO()
    make_send_script([], "sbatch --ntasks [[threads]]", buf)
    content = buf.getvalue()
    assert "[[threads]]" not in content
    assert "1" in content  # default threads value


def test_create_bash_files_nonzero_returncode_logs_error(tmp_path, monkeypatch, caplog):
    infile = tmp_path / "mol.in"
    infile.write_text("")
    monkeypatch.setattr(
        "molbench.bash_wrapper.subprocess.run",
        lambda cmd, shell=False, check=False: subprocess.CompletedProcess(
            args=cmd, returncode=127
        ),
    )
    from molbench.bash_wrapper import create_bash_files

    with caplog.at_level("ERROR", logger="molbench"):
        result = create_bash_files([str(infile)], "definitely_not_a_real_command")
    assert result == []
    assert any("Command failed" in rec.message for rec in caplog.records)


def test_create_bash_files_bare_filename_no_directory(tmp_path, monkeypatch):
    # A file path with no directory component (dirname == "") used to crash
    # os.chdir("") with FileNotFoundError.
    infile = tmp_path / "mol.in"
    infile.write_text("")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", _make_fake_run(".sh"))
    from molbench.bash_wrapper import create_bash_files

    result = create_bash_files(["mol.in"], "fakegen")
    assert len(result) == 1


def test_create_bash_files_restores_cwd_on_exception(tmp_path, monkeypatch):
    infile = tmp_path / "mol.in"
    infile.write_text("")

    def raising_run(cmd, shell=False, check=False):
        raise RuntimeError("boom")

    monkeypatch.setattr("molbench.bash_wrapper.subprocess.run", raising_run)
    original_cwd = os.getcwd()
    from molbench.bash_wrapper import create_bash_files

    with pytest.raises(RuntimeError):
        create_bash_files([str(infile)], "fakegen")
    assert os.getcwd() == original_cwd


def test_make_send_script_uses_abs_path(tmp_path):
    from molbench.bash_wrapper import make_send_script

    f = tmp_path / "job.sh"
    buf = io.StringIO()
    make_send_script([str(f)], "sbatch", buf)
    # The folder path should appear in the script
    assert str(tmp_path) in buf.getvalue()
