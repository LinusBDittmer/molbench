"""Configuration main class.

"""

import os
import json
import molbench.logger as log


class Configuration(dict):
    """
    Class representing a configuration dictionary.

    Attributes
    ----------
    required_fields : dict
        Dictionary containing required configuration fields and their default values.

    Methods
    -------
    load_from_file()
        Load configuration from a JSON file specified by environment variable or default path.
    check_required_fields()
        Check if all required fields are present in the configuration and fill in missing ones 
        with default values.
    
    """

    required_fields = {
        "threads": 1,
        "memory": 50000,
        "walltime": "12:00:00"
    }

    def __init__(self, *args, **kwargs):
        """
        Initialize the Configuration object.

        Parameters
        ----------
        *args, **kwargs : 
            Arguments to initialize the superclass dict.

        """
        self.load_from_file()
        super().__init__(self, *args, **kwargs)

    def load_from_file(self):
        """
        Load configuration from a JSON file.

        This method attempts to load the configuration from a JSON file specified by the 
        environment variable "MOLBENCH_CONFIG" or from the default path "local_config.json" 
        relative to the current file.

        """
        current_dir = os.path.dirname(os.path.realpath(__file__))
        default_config = os.path.join(current_dir, "local_config.json")
        config_path = os.environ.get("MOLBENCH_CONFIG", default_config)
        try:
            with open(config_path, "r") as f:
                self.update(json.load(f))
        except Exception:
            log.critical(f"Configuration file at {config_path} could not be "
                         "parsed.", self)

    def check_required_fields(self):
        """
        Check and fill in missing required fields with default values.

        This method iterates through the required fields defined in the class and checks if they 
        are present in the configuration dictionary. If any required field is missing, it is added
        to the configuration with its default value.

        """
        for field, val in self.required_fields.items():
            if field not in self:
                log.warning(f"Expected Configuration value {field} to be set. "
                            f"Reverting to hardcoded standard of {val}",
                            self)
                self[field] = val

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
