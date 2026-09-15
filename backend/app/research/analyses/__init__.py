"""Analysis implementations. Importing this module registers them with the job
queue, so `jobs.registered()` reflects what can actually run.

Keep this list in step with `capability.assess` — an analysis the matrix marks
available but that is not registered here produces a card the user can click
and nothing else.
"""
from . import (association, burden, descriptive, gwas,  # noqa: F401
               prs, survival)
