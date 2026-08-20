"""Configuration main class."""

import json
import os
from typing import ClassVar

import molbench.logger as log


class Configuration(dict):
    """
    Class representing a configuration dictionary.

    Attributes
    ----------
    required_fields : dict
        Dictionary containing required configuration fields and their default
        values.

    Methods
    -------
    load_from_file()
        Load configuration from a JSON file specified by environment variable
        or default path.

    """

    required_fields: ClassVar = {"threads": 1, "memory": 50000, "walltime": "12:00:00"}

    def __init__(self, *args, **kwargs):
        """
        Initialize the Configuration object.

        Parameters
        ----------
        *args, **kwargs :
            Arguments to initialize the superclass dict.

        """
        self.load_from_file()
        super().__init__(*args, **kwargs)

    def load_from_file(self):
        """
        Load configuration from a JSON file.

        This method attempts to load the configuration from a JSON file
        specified by the environment variable "MOLBENCH_CONFIG" or from the
        default path "local_config.json" relative to the current file.

        """
        current_dir = os.path.dirname(os.path.realpath(__file__))
        default_config = os.path.join(current_dir, "local_config.json")
        config_path = os.environ.get("MOLBENCH_CONFIG", default_config)
        try:
            with open(config_path, "r") as f:
                self.update(json.load(f))
        except Exception:  # noqa: BLE001
            log.critical(
                f"Configuration file at {config_path} could not be parsed.",
                "Configuration",
            )

    def __setattr__(self, attr: str, val) -> None:
        """
        Set attribute value.

        This method sets the value of the specified attribute.

        Parameters
        ----------
        attr : str
            Attribute name.
        val :
            Attribute value.

        """
        self[attr] = val


config = Configuration()
