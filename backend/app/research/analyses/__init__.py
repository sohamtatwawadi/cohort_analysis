"""Analysis implementations. Importing this module registers them with the job
queue, so `jobs.registered()` reflects what can actually run."""
from . import association, burden, gwas, prs, survival  # noqa: F401
