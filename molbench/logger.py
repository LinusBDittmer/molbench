import logging
import os
import sys


instance = None


class MolbenchFormatter(logging.Formatter):
    COLORS = {
        'DEBUG': '\033[1;32m',  # Green
        'INFO': '',  # No color
        'WARNING': '\033[3;33m',  # Orange
        'ERROR': '\033[1;31m',  # Red
        'CRITICAL': '\033[1;35m'  # Pink
    }
    RESET_SEQ = '\033[0m'

    def format(self, record):
        levelname = record.levelname
        if levelname in self.COLORS:
            color_prefix = self.COLORS[levelname]
            record.msg = f"{color_prefix}{record.msg}{self.RESET_SEQ}"
        return super().format(record)


def __init_log_instance():
    global instance
    if instance is not None:
        return
    instance = logging.getLogger("molbench")
    instance.setLevel(os.environ.get("MOLBENCH_VERBOSE", "INFO"))
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(os.environ.get("MOLBENCH_VERBOSE", "INFO"))
    formatter = MolbenchFormatter()
    stream_handler.setFormatter(formatter)
    instance.addHandler(stream_handler)


def _format_cause(cause):
    if isinstance(cause, str):
        return cause
    return str(cause)


def debug(msg: str, cause=None):
    fcause = _format_cause(cause)
    global instance
    if cause is None:
        instance.debug(msg)
    else:
        instance.debug(f"[{fcause}] {msg}")


def info(msg: str, cause=None):
    fcause = _format_cause(cause)
    global instance
    if cause is None:
        instance.info(msg)
    else:
        instance.info(f"[{fcause}] {msg}")


def warning(msg: str, cause=None):
    fcause = _format_cause(cause)
    global instance
    if cause is None:
        instance.warning(msg)
    else:
        instance.warning(f"[{fcause}] {msg}\n")


def error(msg: str, cause=None, etype: str = ""):
    fcause = _format_cause(cause)
    global instance
    if cause is None:
        instance.error(f"{msg} (Error type: {etype})")
    else:
        instance.error(f"[{fcause}] {msg} (Error type: {etype})\n")


def critical(msg: str, cause):
    fcause = _format_cause(cause)
    global instance
    instance.critical(f"[{fcause}] CRITICAL ERROR: {msg}\n")
    sys.exit(-1)


__init_log_instance()
