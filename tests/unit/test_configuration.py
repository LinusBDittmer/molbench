import json

import pytest

from molbench.configuration import Configuration


@pytest.fixture
def cfg_file(tmp_path, monkeypatch):
    f = tmp_path / "cfg.json"
    f.write_text(json.dumps({"threads": 1}))
    monkeypatch.setenv("MOLBENCH_CONFIG", str(f))
    return f


def test_configuration_accepts_positional_mapping(cfg_file):
    # super().__init__(self, *args, **kwargs) used to pass `self` twice,
    # crashing dict.__init__ on any positional argument.
    cfg = Configuration({"memory": 4000})
    assert cfg["memory"] == 4000
    assert cfg["threads"] == 1


def test_configuration_kwargs_still_work(cfg_file):
    cfg = Configuration(memory=8000)
    assert cfg["memory"] == 8000


def test_configuration_missing_file_exits(tmp_path, monkeypatch):
    monkeypatch.setenv("MOLBENCH_CONFIG", str(tmp_path / "nonexistent.json"))
    with pytest.raises(SystemExit):
        Configuration()
