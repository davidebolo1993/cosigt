#
# Based on lsf CookieCutter.py
#
import os
import json

d = os.path.dirname(__file__)
with open(os.path.join(d, "settings.json")) as fh:
    settings = json.load(fh)


def from_entry_or_env(values, key):
    """Return value from ``values`` and override with environment variables."""
    if key in os.environ:
        return os.environ[key]
    else:
        return values[key]


class CookieCutter:

    SBATCH_DEFAULTS = from_entry_or_env(settings, "SBATCH_DEFAULTS")
    CLUSTER_NAME = from_entry_or_env(settings, "CLUSTER_NAME")
    CLUSTER_CONFIG = from_entry_or_env(settings, "CLUSTER_CONFIG")

    @staticmethod
    def get_cluster_option() -> str:
        cluster = CookieCutter.CLUSTER_NAME
        if cluster != "":
            return f"--cluster={cluster}"
        return ""

    @staticmethod
    def get_cluster_logpath() -> str:
        return "logs/slurm/%r/%j"

    @staticmethod
    def get_cluster_jobname() -> str:
        return "%r_%w"
