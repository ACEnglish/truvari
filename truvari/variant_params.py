"""
Truvari main parameters
"""
import types


class VariantParams():
    """
    Holds variant parsing and matching parameters.

    Example
        >>> import truvari
        >>> p = truvari.VariantParams(pctseq=0)
        >>> v = truvari.VariantFile('repo_utils/test_files/variants/input1.vcf.gz', params=p)
        >>> one = next(v); two = next(v)
        >>> one.match(two)
        <truvari.bench.MatchResult (False 2.381)>

    Look at `VariantParams.make_params()` for a list of all params and their defaults
    """

    def __init__(self, args=None, **kwargs):
        """
        Initalize. args is a Namespace from argparse
        """
        if args is not None:
            params = self.make_params_from_args(args)
        else:
            params = self.make_params()

        # Override parameters with those provided in kwargs
        for key, value in kwargs.items():
            if hasattr(params, key):
                setattr(params, key, value)
            else:
                raise ValueError(f"Invalid parameter: {key}")

        for key, value in params.__dict__.items():
            setattr(self, key, value)

    @staticmethod
    def make_params():
        """
        Makes a simple namespace of matching parameters. Holds defaults
        """
        params = types.SimpleNamespace()
        params.reference = None
        params.refdist = 500
        params.pctseq = 0.70
        params.pctsize = 0.70
        params.pctovl = 0.0
        params.typeignore = False
        params.no_roll = False
        params.chunksize = 1000
        params.bSample = 0
        params.cSample = 0
        params.dup_to_ins = False
        params.bnddist = 100
        params.sizemin = 50
        params.sizefilt = 30
        params.sizemax = 50000
        params.passonly = False
        params.no_ref = False
        params.pick = 'single'
        params.ignore_monref = True
        params.check_multi = True
        params.check_monref = True
        params.no_single_bnd = True
        params.write_resolved = False
        params.short_circuit = False
        params.skip_gt = False
        return params

    @staticmethod
    def make_params_from_args(args):
        """
        Makes a simple namespace of matching parameters.
        Populates defaults from make_params, then updates with values from args.
        """
        ret = VariantParams.make_params()

        for key in vars(ret):
            if hasattr(args, key):
                setattr(ret, key, getattr(args, key))

        return ret
