"""
Truvari main parameters
"""

class VariantParams():
    """
    Holds variant parsing and matching parameters.

    Attributes
    ----------
    .. list-table::
       :header-rows: 1

       * - Attribute
         - Description
       * - `refdist`
         - Distance threshold for comparing positions in the reference genome. Default: 500.
       * - `pctseq`
         - Minimum percentage of sequence similarity required for a match. Default: 0.70 (70%).
       * - `pctsize`
         - Minimum percentage of size similarity required for a match. Default: 0.70 (70%).
       * - `pctovl`
         - Minimum percentage of reciprocal overlap required for comparing variants. Default: 0.0 (disabled).
       * - `typeignore`
         - Whether to ignore variant type mismatches during comparison. Default: `False`.
       * - `no_roll`
         - Disable rolling of sequences for comparisons. Default: `False`.
       * - `chunksize`
         - Number of entries to process in each chunk. Default: 1000.
       * - `bSample`
         - Sample name or index for the "base" (a.k.a. self) variants during comparisons. Default: 0.
       * - `cSample`
         - Sample name or index for the "comparison" (a.k.a. other) variants during comparisons. Default: 0.
       * - `dup_to_ins`
         - Whether to treat duplications as insertions for some operations. Default: `False`.
       * - `bnddist`
         - Maximum allowed distance for breakend (BND) comparisons. Default: 100.
       * - `sizemin`
         - Minimum variant size to consider. Default: 50.
       * - `sizefilt`
         - Minimum size filter for comparison in the "comparison" dataset. Default: 30.
       * - `sizemax`
         - Maximum variant size to consider. Default: 50,000; -1, off.
       * - `passonly`
         - Whether to only consider variants with a "PASS" filter status. Default: `False`.
       * - `no_ref`
         - Whether to ignore reference homozygous variants in (a)ll, (b)ase, or (c)omp VCF Default: `False` (off)`.
       * - `pick`
         - Strategy for picking matches by Bench (single, ac, multi).
       * - `ignore_monref`
         - Whether to ignore monoallelic reference calls. Default: `True`.
       * - `check_multi`
         - Whether to check for and handle multi-allelic records. Default: `True`.
       * - `check_monref`
         - Whether to check for monoallelic reference calls. Default: `True`.
       * - `no_single_bnd`
         - Whether to exclude single-end breakends (BNDs) from comparisons. Default: `True`.
       * - `write_resolved`
         - Whether to write resolved REF/ALT sequences to output. Default: `False`.
       * - `decompose`
         - Decompose symbolic variants (e.g. <DEL>) into BNDs when comparing with BNDs. Default: `True`
       * - `short_circuit`
         - Whether to enable short-circuit logic for early exits in comparisons. Default: `False`.
       * - `skip_gt`
         - Whether to skip genotype comparisons. Default: `False`.
       * - `max_resolve`
         - Maximum size of a variant to attempt sequence resolving. Default: 25000.
       * - `dynthresh`
         - Tuple of parameters for dynamic thresholding (min_diff, max_diff, sz_min, sz_max) Default: `None`=off
    """

    DEFAULTS = {
        "reference": None,
        "refdist": 500,
        "pctseq": 0.70,
        "pctsize": 0.70,
        "pctovl": 0.0,
        "typeignore": False,
        "no_roll": False,
        "chunksize": 1000,
        "bSample": 0,
        "cSample": 0,
        "dup_to_ins": False,
        "bnddist": 100,
        "sizemin": 50,
        "sizefilt": 30,
        "sizemax": 50000,
        "passonly": False,
        "no_ref": False,
        "pick": "single",
        "ignore_monref": True,
        "check_multi": True,
        "check_monref": True,
        "no_single_bnd": True,
        "write_resolved": False,
        "decompose": True,
        "short_circuit": False,
        "skip_gt": False,
        "max_resolve": 25000,
        "dynthresh": None
    }

    def __init__(self, args=None, **kwargs):
        """
        Initialize VariantParams with defaults, args, and kwargs.

        Parameters
        ----------
        args : Namespace (optional)
            An argparse.Namespace object to initialize parameters.
        kwargs : dict
            Additional parameters to override defaults.
        """
        # Start with defaults
        params = self.DEFAULTS.copy()

        # Override with args if provided
        if args:
            for key in vars(args):
                if key in params:
                    params[key] = getattr(args, key)

        # Override with kwargs
        for key, value in kwargs.items():
            if key in params:
                params[key] = value
            else:
                raise ValueError(f"Invalid parameter: {key}")

        # Set attributes
        self.__dict__.update(params)

    def calc_dyn_thresh(self, size):
        """
        # min_diff=5, max_diff=30, s_min=50, s_max=1500
        Linear interpolation of size to determine threshold
        """
        min_diff, max_diff, s_min, s_max = self.dynthresh
        if s_min is None:
            s_min = 2 * min_diff
        if s_max is None:
            s_max = 2 * max_diff

        if size <= s_min:
            allowed_diff = min_diff
        elif size >= s_max:
            allowed_diff = max_diff
        else:
            # Linear interpolation between min_diff and max_diff
            scale = (size - s_min) / (s_max - s_min)
            allowed_diff = min_diff + scale * (max_diff - min_diff)

        return 1 - allowed_diff / size

    def get_pctsize(self, size):
        """
        Returns either the flat self.pctsize or the dynamic threshold based on the point
        """
        return self.calc_dyn_thresh(size) if self.dynthresh else self.pctsize

    def get_pctseq(self, size):
        """
        Returns either the flat self.pctseq or the dynamic threshold based on the point
        """
        return self.calc_dyn_thresh(size) if self.dynthresh else self.pctseq
