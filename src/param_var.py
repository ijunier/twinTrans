import argparse
import os, sys
import numpy as np
import shutil

"""
    Parameters and variables used by ../bin/twin.py
    For lightness, only downstream quantities are referred to as _down
    i.e., upstream quantities are NOT referred to as _up -- this might change in the future
"""


def parsing_cmd():

    parser = argparse.ArgumentParser(
        description="Run the twin transcriptional-loop model in presence of topoisomerases",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # POSITIONAL ARGUMENTS
    parser.add_argument(
        "output_folder", help="results output folder (-f to force overwriting)"
    )

    # OPTIONAL ARGUMENTS
    parser.add_argument(
        "-f", action="store_true", help="overwrite results output folder\n"
    )
    parser.add_argument(
        "-promfollow", action="store_true", help="following promoter status\n"
    )

    # GENE
    parser.add_argument(
        "-tss",
        metavar="X",
        dest="gene_tss",
        type=int,
        help="distance between TSS and upstream barrier (in bp)",
        default=55,
    )
    parser.add_argument(
        "-L",
        metavar="X",
        dest="gene_L",
        type=int,
        help="gene length (in bp)",
        default=900,
    )
    parser.add_argument(
        "-Ld",
        metavar="X",
        dest="gene_Ldown",
        type=int,
        help="distance seperating gene termination STOP and downstream barrier (in bp)",
        default=320,
    )

    # PROMOTER
    parser.add_argument(
        "-kb",
        metavar="X",
        dest="promoter_kb",
        type=float,
        help="promoter binding rate (per s)",
        default=1,
    )
    parser.add_argument(
        "-ko",
        metavar="X",
        dest="promoter_ko",
        type=float,
        help="OC formation rate (per s)",
        default=1,
    )
    parser.add_argument(
        "-so",
        metavar="X",
        dest="promoter_sigma_o",
        type=float,
        help="supercoiling density threshold for OC formation (no unit)",
        default=-0.05,
    )
    parser.add_argument(
        "-ke",
        metavar="X",
        dest="promoter_ke",
        type=float,
        help="OC escape rate (per s)",
        default=1,
    )

    # RNAP
    parser.add_argument(
        "-vm",
        metavar="X",
        dest="rnap_vm",
        type=float,
        help="RNAP maximum translocation speed (in bp per s)",
        default=25,
    )
    parser.add_argument(
        "-G",
        metavar="X",
        dest="rnap_torque_stall",
        type=float,
        help="RNAP stalling torque (in pN.nm)",
        default=18.5,
    )
    parser.add_argument(
        "-rel",
        metavar="X",
        dest="rnap_excluded_length",
        type=float,
        help="RNAP excluded length (in bp)",
        default=30,
    )

    # TOPOI
    parser.add_argument(
        "-lansT",
        metavar="X",
        dest="topoI_lambda_ns",
        type=float,
        help="non-specific rate of upstream TopoI activity (per bp per s)",
        default=1e-4,
    )
    parser.add_argument(
        "-LasT",
        metavar="X",
        dest="topoI_Lambda_s",
        type=float,
        help="specific rate of upstream TopoI activity (per s)",
        default=1,
    )
    parser.add_argument(
        "-sa",
        metavar="X",
        dest="topoI_sigma_active",
        type=float,
        help="sigma below which TopoI is active (no unit)",
        default=-0.05,
    )

    # GYRASE
    parser.add_argument(
        "-lansG",
        metavar="X",
        dest="gyrase_lambda_ns",
        type=float,
        help="non-specific rate of downstream gyrase activity (per bp per s)",
        default=1e-4,
    )
    parser.add_argument(
        "-LasG",
        metavar="X",
        dest="gyrase_Lambda_s",
        type=float,
        help="specific rate (absolute value) of downstream gyrase activity (per s)",
        default=2,
    )

    # DNA
    parser.add_argument(
        "-A",
        metavar="X",
        dest="dna_A",
        type=float,
        help="proportionality constant between torque and sigma for unstructured DNA (in pN.nm)",
        default=300,
    )

    # SIMULATION
    parser.add_argument(
        "-Ni",
        metavar="X",
        dest="simup_Niterations",
        type=int,
        help="maximum number of iterations",
        default=int(1e10),
    )
    parser.add_argument(
        "-Nt",
        metavar="X",
        dest="simup_Ntranscripts_max",
        type=int,
        help="maximum number of transcripts",
        default=int(1e4),
    )
    parser.add_argument(
        "-Net",
        metavar="X",
        dest="simup_Nevery_transcripts",
        type=int,
        help="writing out results every X",
        default=int(1e3),
    )

    # RESOLUTION
    parser.add_argument(
        "-dx",
        metavar="X",
        dest="coarse_g_dx",
        type=float,
        help="resolution, i.e., RNAP translocation step (in bp)",
        default=1,
    )

    # DYNAMICS
    parser.add_argument(
        "-dyn",
        choices=["up2down", "async", "down2up"],
        dest="model_dyn",
        help="types of dynamics for the stochastic updating",
        default="async",
    )

    # TESTING/DEBUGGING
    parser.add_argument(
        "-verbose",
        "-v",
        action="store_true",
        help="(for testing) provides some output of the numbers",
    )
    parser.add_argument(
        "-v_every",
        metavar="X",
        dest="simup_verbose_every",
        type=int,
        help="verbose every X iterations (-verbose is required)",
        default=int(1e3),
    )
    parser.add_argument(
        "-rseed",
        metavar="X",
        type=int,
        help="(for debugging) specifies the seed for random number generator",
    )

    args = parser.parse_args()

    if args.rseed:
        np.random.seed(seed=args.rseed)
    # fixing random seed

    # OUTPUT FOLDER
    if args.f and os.path.isdir(args.output_folder):
        shutil.rmtree(args.output_folder, ignore_errors=True)
    try:
        os.mkdir(args.output_folder)
    except OSError as err:
        print("\n" + str(err))
        sys.exit("(use -f to overwrite it)\n")

    return args


class SimuParam:
    """Simulation parameters"""

    def __init__(self, args):

        self.verbose = args.verbose
        self.verbose_every = args.simup_verbose_every
        self.promoter_to_follow = args.promfollow

        self.Niterations = args.simup_Niterations
        # total number of iterations
        self.Ntranscripts_max = args.simup_Ntranscripts_max
        # maximum number of transcripts
        self.Nevery_transcripts = args.simup_Nevery_transcripts
        # writing out every Nevery_transcripts

        self.fo_out = args.output_folder  # output folder


class ModelParam:
    """Modelling parameters"""

    def __init__(self, args):

        self.t_dyn = args.model_dyn
        # dynamics types
        # asynchronous: async (default)
        # up2down: up2down (from upstream RNAP to downstream RNAP)
        # down2up: down2up (from downstream RNAP to upstream RNAP)

        self.dna = self._DNA(args)
        self.gene = self._Gene(args, self.dna)
        self.promoter = self._Promoter(args)
        self.rnap = self._RNAP(args, self.dna)
        self.topoI = self._TopoI(args)
        self.gyrase = self._Gyrase(args)
        self.coarse_g = self._CoarseGraining(
            args, self.dna, self.rnap, self.topoI, self.gyrase
        )

    class _DNA:
        """DNA microscopic parameters"""

        def __init__(self, args):

            self.n = 10.5  # bp
            # structural properties

            self.A = args.dna_A  # pN.nm
            # torque <-> sigma properties

    class _Gene:
        """Gene parameters, exlcuding promoter features"""

        def __init__(self, args, dna):

            self.tss = args.gene_tss
            # TSS position
            self.rnap_xi = self.tss - int(args.rnap_excluded_length / 2)
            # initial position of rnap at the promoter
            self.Lc_lk_rnap_xi = self.rnap_xi / dna.n
            # TSS distance to upstream barrier in units of linking humber

            self.L = args.gene_L
            # gene length (in bp)
            self.Ldown = args.gene_Ldown
            # downstream distance to barrier (in bp)

            self.Ldomain = self.tss + self.L + self.Ldown
            # total domain length
            self.Ldomain_lk_relaxed = self.Ldomain / dna.n
            # in unit of linking number in the relaxed state
            self.Ldomain_lk = self.Ldomain_lk_relaxed
            # current state

    class _Promoter:
        """Promoter features"""

        def __init__(self, args):

            # BINDING RATE
            self.kb = args.promoter_kb
            self.kb_s = 1 / self.kb
            # corresponding scale (s)

            # FORMATION OF THE OPEN COMPLEX
            self.ko = args.promoter_ko
            self.ko_s = 1 / self.ko
            # corresponding scale (s)
            self.sigma_o = args.promoter_sigma_o
            # supercoiling density threshold

            # OC ESCAPE RATE
            self.ke = args.promoter_ke
            self.ke_s = 1 / self.ke
            # corresponding scale (s)

    class _RNAP:
        """RNAP microscopic parameters"""

        def __init__(self, args, dna):

            self.excluded_length_elongation = args.rnap_excluded_length
            # in bp

            self.vm = args.rnap_vm
            # max velocity of the RNAP (bp per s)
            self.torque_stall = args.rnap_torque_stall
            # 18.5 pN in the presence of GreB (in vivo situation) (pN.nm)

            self.sigma_stall = -self.torque_stall / dna.A
            # (sigma's are signed here)

    class _TopoI:
        """TopoI activity"""

        def __init__(self, args):

            self.sigma_active = args.topoI_sigma_active
            # sigma below which TopoI is active

            self.lambda_ns = args.topoI_lambda_ns
            # non-specific (per base pair and per second)

            self.Lambda_s = args.topoI_Lambda_s
            # specific (~ at the promoter)

    class _Gyrase:
        """Gyrase activity
        rem: it is active whenever sigma is larger than stalling sigma's (i.e., stalling torque of the RNAP)
        """

        def __init__(self, args):

            self.lambda_ns = args.gyrase_lambda_ns
            # non-specific (per base pair and per second)

            self.Lambda_s = args.gyrase_Lambda_s
            # specific (~ at the most downstream RNAP)

    class _CoarseGraining:
        """Coarse-graining parameters (modulated by resolution, i.e. dx)"""

        def __init__(self, args, dna, rnap, topoI, gyrase):

            self.dx = args.coarse_g_dx
            # translocation step (in bp)
            self.dLk = self.dx / dna.n
            # translocation step in Lk unit
            self.tau_0 = self.dx / rnap.vm
            # corresponding time unit

            # TOPOI
            self.p_topoI_ns_per_bp = topoI.lambda_ns * self.tau_0
            # maximal probability to non-specifically generate +1 Lk per bp
            self.p_topoI_s = np.min((topoI.Lambda_s * self.tau_0, 1))
            # probability to specifically generate it (~ at the promoter)

            # GYRASE
            self.p_gyrase_ns_per_bp = gyrase.lambda_ns * self.tau_0 / 2
            # maximal probability to non-specifically generate -2 Lk per bp
            self.p_gyrase_s = np.min((gyrase.Lambda_s * self.tau_0 / 2, 1))
            # probability to specifically generate it (~ at the downstream RNAP)

    def _test(self):
        """Some tests"""

        # RNAP
        if self.rnap.sigma_stall > 0:
            sys.exit("stalling sigma's must be none-positive")

        # TOPOI
        if self.topoI.sigma_active <= self.rnap.sigma_stall:
            sys.exit("sigma_active <= sigma_stall is not biologically relevant")


class RNAP:
    def __init__(self):

        self.tb = 0
        # time at binding event
        self.tocf = 0
        # time at OC formation event
        self.tesc = 0
        # time at OC escape

        self.t_elongating = False
        # = True if RNAP is elongating

        self.X = 0
        # position in base pairs

        # UPSTREAM QUANTITIES
        # we drop "up" for lightness
        self.Lc_lk = None
        # upstream contour length in units of Lk
        self.Lk = None
        # upstream linking number
        self.sigma = None
        # corresponding supercoiling density
        self.pstruct = None
        # corresponding supercoiling density

        # DOWNSTREAM QUANTITIES
        self.Lc_down_lk = None
        # downstream contour length in units of Lk
        self.Lk_down = None
        # downstream linking number
        self.sigma_down = None
        # corresponding supercoiling density


class Trajectory:
    def __init__(self):

        self.niter = 0
        # number of iteration
        self.time = 0
        # real time

        self.Ntranscripts = 0
        # number of produced transcripts

        # STATISTICAL ANALYSIS OF TIMES
        self.prod_times = {"n": 0, "mean": 0}
        # time between the production of 2 RNA transcripts
        self.time_last_prod = 0
        self.b2b_times = {"n": 0, "mean": 0}
        # time between two binding events
        self.ref_binding_time = 0
        self.ocf_times = {"n": 0, "mean": 0}
        # time for OC to form once bound
        self.esc_times = {"n": 0, "mean": 0}
        # time to escape promoter once OC is formed
        self.initiation_times = {"n": 0, "mean": 0}
        # total time for initiation once it is bound
        self.elongation_times = {"n": 0, "mean": 0}
        # time for elongation once it is escaped from promoter

        self.sigma_at_ocf = []
        # supercoiling density at of OC formation stage
