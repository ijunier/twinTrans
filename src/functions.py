import numpy as np, os, sys
import math
from param_var import *

"""
    Functions used by ../bin/twin.py
    For lightness, only downstream quantities are referred to as _down
    i.e., upstream quantities are NOT referred to as _up -- this might change in the future
"""

# supercoiling densities
def _sigma(rnap, modelP):
    return (rnap.Lk - rnap.Lc_lk) / rnap.Lc_lk


def _sigma_down(rnap, modelP):
    return (rnap.Lk_down - rnap.Lc_down_lk) / rnap.Lc_down_lk


# Simulation runs
def generate_run_follow_promoter(modelP, simuP):

    traj = Trajectory()
    RNAP_list = []
    # DNA-bound RNAPs

    nextevent2iter = {tag: -1 for tag in ["b", "oc", "esc"]}
    # next trial for binding, OC formation and promoter escape
    # -1 to avoid initial True
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = int(Z / modelP.coarse_g.tau_0)

    write_follow_promoter(traj, None, modelP, simuP, header=True)
    while (
        traj.niter < simuP.Niterations
        and traj.Ntranscripts < simuP.Ntranscripts_max
    ):
        write_follow_promoter(traj, RNAP_list, modelP, simuP)

        # Binding
        if traj.niter == nextevent2iter["b"]:
            binding_stage(modelP, RNAP_list, traj, nextevent2iter)

        # OC formation
        if traj.niter == nextevent2iter["oc"]:
            oc_formation_stage(modelP, RNAP_list, traj, nextevent2iter)

        # Promoter escape
        if traj.niter == nextevent2iter["esc"]:
            escape_stage(modelP, RNAP_list, traj)

        # Topoisomerases
        topo_stage(modelP, RNAP_list)

        # Elongation
        elongation_stage(modelP, RNAP_list)

        # Termination
        termination_stage(modelP, simuP, RNAP_list, traj)

        traj.niter += 1
        traj.time = traj.niter * modelP.coarse_g.tau_0

    return


def generate_run_multiple_transcrtipts(modelP, simuP):

    traj = Trajectory()
    RNAP_list = []
    # DNA-bound RNAPs

    nextevent2iter = {tag: -1 for tag in ["b", "oc", "esc"]}
    # nextevent2iter: next trial for binding, OC formation and promoter escape
    # -1 to avoid initial True
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = int(Z / modelP.coarse_g.tau_0)

    write_transcripts_on_the_fly(traj, simuP, header=True)
    while (
        traj.niter < simuP.Niterations
        and traj.Ntranscripts < simuP.Ntranscripts_max
    ):
        verbosing(simuP, traj, RNAP_list)

        # Binding
        if traj.niter == nextevent2iter["b"]:
            binding_stage(modelP, RNAP_list, traj, nextevent2iter)

        # OC formation
        if traj.niter == nextevent2iter["oc"]:
            oc_formation_stage(modelP, RNAP_list, traj, nextevent2iter)

        # Promoter escape
        if traj.niter == nextevent2iter["esc"]:
            escape_stage(modelP, RNAP_list, traj)

        # Topoisomerases
        topo_stage(modelP, RNAP_list)

        # Elongation
        elongation_stage(modelP, RNAP_list)

        # Termination
        termination_stage(modelP, simuP, RNAP_list, traj)

        traj.niter += 1
        traj.time = traj.niter * modelP.coarse_g.tau_0

    return traj


# Transcription stages
def binding_stage(modelP, RNAP_list, traj, nextevent2iter):
    """Binding of the RNAP at the promoter"""

    if (
        not RNAP_list
        or RNAP_list[-1].X > modelP.gene.rnap_xi + modelP.rnap.excluded_length_elongation
    ):
        # promoter is free => binding occurs!

        if traj.ref_binding_time > 0:
            traj.b2b_times["mean"] = (
                traj.b2b_times["n"] * traj.b2b_times["mean"]
                + traj.time
                - traj.ref_binding_time
            ) / (traj.b2b_times["n"] + 1)
            traj.b2b_times["n"] += 1

        traj.ref_binding_time = traj.time

        # SETTING UP NEW RNAP
        rnap = RNAP()

        rnap.tb = traj.time
        rnap.t_elongating = False
        rnap.X = modelP.gene.rnap_xi
        rnap.Lc_lk = modelP.gene.Lc_lk_rnap_xi

        # SIGMA PROPERTIES                
        if not RNAP_list:
            # a single RNAP => topological properties dictated by the entire topological domain
            rnap.sigma = (
                modelP.gene.Ldomain_lk - modelP.gene.Ldomain_lk_relaxed
            ) / modelP.gene.Ldomain_lk_relaxed
            rnap.Lc_down_lk = modelP.gene.Ldomain_lk_relaxed - modelP.gene.Lc_lk_rnap_xi
        else:
            # multiple RNAP => topological properties dictated by the immediate downstream RNAP
            rnap.sigma = RNAP_list[-1].sigma
            rnap.Lc_down_lk = RNAP_list[-1].Lc_lk - rnap.Lc_lk

        rnap.Lk = (1 + rnap.sigma) * rnap.Lc_lk

        rnap.sigma_down = rnap.sigma
        rnap.Lk_down = (1 + rnap.sigma_down) * rnap.Lc_down_lk

        # ADDING THE NEW RNAP
        RNAP_list.append(rnap)

        # NEXT STAGE: OC FORMATION
        Z = np.random.exponential(scale=modelP.promoter.ko_s, size=None)
        nextevent2iter["oc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1

    # NEXT BINDING TRIAL (INDEPENDENT WHETHER BINDING HAS OCCURED OR NOT)
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = np.max(
        (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
    )  # at least niter + 1

    return


def oc_formation_stage(modelP, RNAP_list, traj, nextevent2iter):
    """OC formation if sigma <= threshold"""

    traj.sigma_at_ocf.append(RNAP_list[-1].sigma)
    if RNAP_list[-1].sigma <= modelP.promoter.sigma_o:
        # sigma is below threshold: OC formation occurs!

        RNAP_list[-1].tocf = traj.time

        traj.ocf_times["mean"] = (
            traj.ocf_times["n"] * traj.ocf_times["mean"]
            + RNAP_list[-1].tocf
            - RNAP_list[-1].tb
        ) / (traj.ocf_times["n"] + 1)
        traj.ocf_times["n"] += 1

        # NEXT STAGE: ESCAPE EVENT
        Z = np.random.exponential(scale=modelP.promoter.ke_s, size=None)
        nextevent2iter["esc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1
    else:
        # sigma is above threshold: OC formation has failed
        
        # NEXT OC FORMATION TRIAL
        Z = np.random.exponential(scale=modelP.promoter.ko_s, size=None)
        nextevent2iter["oc"] = np.max(
            (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
        )  # at least niter + 1

    return


def escape_stage(modelP, RNAP_list, traj):
    """promoter escape => RNAP is now in elongating mode"""

    RNAP_list[-1].t_elongating = True
    RNAP_list[-1].tesc = traj.time

    traj.esc_times["mean"] = (
        traj.esc_times["n"] * traj.esc_times["mean"]
        + RNAP_list[-1].tesc
        - RNAP_list[-1].tocf
    ) / (traj.esc_times["n"] + 1)
    traj.esc_times["n"] += 1

    traj.initiation_times["mean"] = (
        traj.initiation_times["n"] * traj.initiation_times["mean"]
        + RNAP_list[-1].tesc
        - RNAP_list[-1].tb
    ) / (traj.initiation_times["n"] + 1)
    traj.initiation_times["n"] += 1

    # UPDATING DOWNSTREAM RNAP
    # change Lk because the new elongating RNA becomes a barrier
    # one can actually check the conservation of the Lk
    if len(RNAP_list) > 1:
        RNAP_list[-2].Lc_lk -= RNAP_list[-1].Lc_lk
        RNAP_list[-2].Lk = (1 + RNAP_list[-2].sigma) * RNAP_list[-2].Lc_lk

    return


def topo_stage(modelP, RNAP_list):
    """TopoI and gyrase activity"""
    
    if RNAP_list:
        topo_stage_RNAPpresent(modelP, RNAP_list)
    else:
        topo_stage_RNAPabsent(modelP)

    return


def topo_stage_RNAPpresent(modelP, RNAP_list):
    """TopoI and gyrase activity in the presence of at least one DNA-bound RNAP"""

    # UPSTREAM
    # non-specific activities
    DTopoI, DGyrase = 0, 0
    if not RNAP_list[-1].t_elongating:
        # the most upstream RNAP (at the promoter) is not a barrier
        if len(RNAP_list) == 1:
            # topoisomerases can act anywhere along the domain
            domain_length_topo = modelP.gene.Ldomain
            DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
            DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
        else:
            # the second RNAP is a barrier and we consider activity upstream
            domain_length_topo = RNAP_list[-2].Lc_lk * modelP.dna.n
            DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
            # RNAP_list[-1] because RNAP_list[-1].sigma = RNAP_list[-2].sigma here
            DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
            # RNAP_list[-1] because RNAP_list[-1].sigma = RNAP_list[-2].sigma here
    else:
        # the most upstream RNAP is a barrier and we consider activity upstream
        domain_length_topo = RNAP_list[-1].Lc_lk * modelP.dna.n
        DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
        DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)

    # specific activity in the presence of transcription (only TopoI)
    DTopoI_spec = 0
    if len(RNAP_list) > 1 or RNAP_list[-1].t_elongating:
        DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)

    if DTopoI != 0 or DGyrase != 0 or DTopoI_spec != 0:  # updating topo properties
        modelP.gene.Ldomain_lk += DTopoI + DGyrase + DTopoI_spec

        if not RNAP_list[-1].t_elongating:
            # properties of non-elongating RNAP are dictated by its downstream RNAP (if it exists)
            if len(RNAP_list) == 1:
                RNAP_list[0].sigma = (
                    modelP.gene.Ldomain_lk - modelP.gene.Ldomain_lk_relaxed
                ) / modelP.gene.Ldomain_lk_relaxed
            else:
                RNAP_list[-2].Lk += DTopoI + DGyrase + DTopoI_spec
                RNAP_list[-2].sigma = _sigma(RNAP_list[-2], modelP)
                RNAP_list[-1].sigma = RNAP_list[-2].sigma

            RNAP_list[-1].Lk = (1 + RNAP_list[-1].sigma) * RNAP_list[-1].Lc_lk
            RNAP_list[-1].sigma_down = RNAP_list[
                -1
            ].sigma  # because RNAP is not a barrier
            RNAP_list[-1].Lk_down = (1 + RNAP_list[-1].sigma_down) * RNAP_list[
                -1
            ].Lc_down_lk
        else:
            # RNAP is a barrier
            RNAP_list[-1].Lk += DTopoI + DGyrase + DTopoI_spec
            RNAP_list[-1].sigma = _sigma(RNAP_list[-1], modelP)

    # DOWNSTREAM
    DTopoI_down, DGyrase_down, DGyrase_spec = 0, 0, 0
    if RNAP_list[0].t_elongating:
        # if non elongating, this means a single non-elongating RNAP => treated at the upstream level
        domain_length_topo = RNAP_list[0].Lc_down_lk * modelP.dna.n
        DTopoI_down = DLk_TopoI(domain_length_topo, RNAP_list[0], modelP, dna="down")
        DGyrase_down = DLk_Gyrase(domain_length_topo, RNAP_list[0], modelP, dna="down")
        DGyrase_spec = DLk_Gyrase("spec", RNAP_list[0], modelP, dna="down")

        modelP.gene.Ldomain_lk += DTopoI_down + DGyrase_down + DGyrase_spec
        RNAP_list[0].Lk_down += DTopoI_down + DGyrase_down + DGyrase_spec
        RNAP_list[0].sigma_down = _sigma_down(RNAP_list[0], modelP)

    return


def topo_stage_RNAPabsent(modelP):
    """TopoI and gyrase activity in the absence of RNAP"""

    # TopoI
    sigma = (
        modelP.gene.Ldomain_lk - modelP.gene.Ldomain_lk_relaxed
    ) / modelP.gene.Ldomain_lk_relaxed
    modelP.gene.Ldomain_lk += DLk_TopoI_noRNAP(modelP.gene.Ldomain, modelP, sigma)

    # gyrase
    sigma = (
        modelP.gene.Ldomain_lk - modelP.gene.Ldomain_lk_relaxed
    ) / modelP.gene.Ldomain_lk_relaxed
    modelP.gene.Ldomain_lk += DLk_Gyrase_noRNAP(modelP.gene.Ldomain, modelP, sigma)

    return

# Elementary generations of linking numbers
def DLk_TopoI(domain_length_topo, rnap, modelP, dna="up"):
    """
        TopoI activity associated with an RNAP
        - worked for both upstream and downstream the "RNAP convoy"
        - not active if sigma > sigma_active
    """

    if dna == "up":
        sigma = rnap.sigma
    else:
        sigma = rnap.sigma_down

    if sigma > modelP.topoI.sigma_active:
        return 0
    else:
        if domain_length_topo != "spec":
            return np.random.poisson(
                modelP.coarse_g.p_topoI_ns_per_bp * domain_length_topo
            )
        else:
            # idpt of distance
            return np.random.uniform() <= modelP.coarse_g.p_topoI_s


def DLk_Gyrase(domain_length_topo, rnap, modelP, dna="up"):
    """
        Gyrase activity associated with an RNAP
        - worked for both upstream and downstream the "RNAP convoy"
        - not active if sigma < sigma_stall
    """

    if dna == "up":
        sigma = rnap.sigma
    else:
        sigma = rnap.sigma_down

    if sigma < modelP.rnap.sigma_stall:
        # RNAP stalling torque is the gyrase threshold
        return 0
    else:
        if domain_length_topo != "spec":
            return -2 * np.random.poisson(
                modelP.coarse_g.p_gyrase_ns_per_bp * domain_length_topo
            )
        else:
            # idpt of distance
            return -2 * (np.random.uniform() <= modelP.coarse_g.p_gyrase_s)


def DLk_TopoI_noRNAP(domain_length_topo, modelP, sigma):
    """
        TopoI activity in the absence of any RNAP => sigma is specified as an argument
        - not active if sigma > sigma_active
    """

    if sigma > modelP.topoI.sigma_active:
        return 0
    else:
        return np.random.poisson(modelP.coarse_g.p_topoI_ns_per_bp * domain_length_topo)


def DLk_Gyrase_noRNAP(domain_length_topo, modelP, sigma):
    """
        Gyrase activity in the absence of any RNAP => sigma is specified as an argument
        - not active if sigma > sigma_active
    """

    if sigma < modelP.rnap.sigma_stall:
        # stalling torque is the gyrase threshold
        return 0
    else:
        return -2 * np.random.poisson(
            modelP.coarse_g.p_gyrase_ns_per_bp * domain_length_topo
        )


def elongation_stage(modelP, RNAP_list):
    """RNAPs translocation"""

    if RNAP_list:
        ix_ = np.arange(len(RNAP_list))

        if modelP.t_dyn == "async":
            np.random.shuffle(ix_)
        elif modelP.t_dyn == "up2down":
            ix_ = ix_[::-1]
        elif modelP.t_dyn != "down2up":
            sys.exit("Error with modelP.t_dyn, which currently is = %s" % modelP.t_dyn)

        for ix_rnap in ix_:
            RNAP_translocation(ix_rnap, RNAP_list, modelP)

    return


# Translocations and their topological consequences
def RNAP_translocation(ix_rnap, RNAP_list, modelP):
    
    if (
        not RNAP_list[ix_rnap].t_elongating
        or not RNAP_list[ix_rnap].sigma >= modelP.rnap.sigma_stall
        or not RNAP_list[ix_rnap].sigma_down <= np.abs(modelP.rnap.sigma_stall)
    ):
        # not elongating or beyond sigma_stall = no translocation
        return

    # Rem1: linking numbers do not change, supercoiling density does
    # Rem2: no constrain on one RNAP overtaking another one (supercoiling constraints do the job)

    RNAP_list[ix_rnap].X += modelP.coarse_g.dx
    RNAP_list[ix_rnap].Lc_lk += modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma = _sigma(RNAP_list[ix_rnap], modelP)

    RNAP_list[ix_rnap].Lc_down_lk -= modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma_down = _sigma_down(RNAP_list[ix_rnap], modelP)

    # UPDATING DOWNSTREAM RNAP
    if ix_rnap > 0:        
        RNAP_list[ix_rnap - 1].Lc_lk -= modelP.coarse_g.dLk
        RNAP_list[ix_rnap - 1].sigma = _sigma(RNAP_list[ix_rnap - 1], modelP)

    # UPDATING UPSTREAM RNAP
    if ix_rnap < len(RNAP_list) - 1:
        if RNAP_list[ix_rnap + 1].t_elongating:
            # the RNAP is elongating        
            RNAP_list[ix_rnap + 1].Lc_down_lk += modelP.coarse_g.dLk
            RNAP_list[ix_rnap + 1].sigma_down = _sigma_down(RNAP_list[ix_rnap + 1], modelP)
        elif ix_rnap == len(RNAP_list) - 2 and not RNAP_list[len(RNAP_list) - 1].t_elongating:
            # the RNAP is NON-elongating
            RNAP_list[-1].sigma = RNAP_list[-2].sigma
            RNAP_list[-1].Lk = (1 + RNAP_list[-1].sigma) * RNAP_list[-1].Lc_lk

            RNAP_list[-1].Lc_down_lk += modelP.coarse_g.dLk
            RNAP_list[-1].sigma_down = RNAP_list[-2].sigma
            RNAP_list[-1].Lk_down = (1 + RNAP_list[-1].sigma_down) * RNAP_list[
                -1
            ].Lc_down_lk

    return

def termination_stage(modelP, simuP, RNAP_list, traj):
    """Termination stage: transcript production by the most downstream RNAP"""

    if RNAP_list and RNAP_list[0].X >= modelP.gene.tss + modelP.gene.L:
        traj.Ntranscripts += 1

        # ELONGATION TIME
        traj.elongation_times["mean"] = (
            traj.elongation_times["n"] * traj.elongation_times["mean"]
            + traj.time
            - RNAP_list[0].tesc
        ) / (traj.elongation_times["n"] + 1)
        traj.elongation_times["n"] += 1

        # PRODUCTION TIME
        if traj.time_last_prod > 0:
            traj.prod_times["mean"] = (
                traj.prod_times["n"] * traj.prod_times["mean"]
                + traj.time
                - traj.time_last_prod
            ) / (traj.prod_times["n"] + 1)
            traj.prod_times["n"] += 1
        traj.time_last_prod = traj.time

        # TRANSCRIPT TERMINATION => RNAP REMOVAL
        del RNAP_list[0]

        if not traj.Ntranscripts % simuP.Nevery_transcripts:
            write_transcripts_on_the_fly(traj, simuP)

    return


# I/O
def verbosing(simuP, traj, RNAP_list):
    """some info to STDOUT (verbosing mode, e.g., for debugging)"""

    if simuP.verbose and not traj.niter % simuP.verbose_every:
        if not len(RNAP_list):
            print(traj.niter, traj.Ntranscripts, len(RNAP_list), end="\r")
        else:
            if len(RNAP_list) == 1:
                print(
                    traj.niter,
                    traj.Ntranscripts,
                    len(RNAP_list),
                    RNAP_list[-1].sigma,
                    RNAP_list[-1].sigma_down,
                    end="\r",
                )
            else:
                print(
                    traj.niter,
                    traj.Ntranscripts,
                    len(RNAP_list),
                    RNAP_list[-1].sigma,
                    RNAP_list[-1].sigma_down,
                    RNAP_list[-1].X,
                    RNAP_list[-1].t_elongating,
                    RNAP_list[-2].sigma,
                    RNAP_list[-2].sigma_down,
                    RNAP_list[-2].X,
                    RNAP_list[-2].t_elongating,
                    RNAP_list[0].sigma_down,
                    RNAP_list[0].X,
                    end="\r",
                )


def write_follow_promoter(traj, RNAP_list, modelP, simuP, header=False):
    """promoter properties"""

    fi = simuP.fo_out + "/traj_promoter.txt"

    if header:
        with open(fi, "w") as out:
            out.write("time\tsigma\tt_upRNAPelong\n")
        return

    with open(fi, "a") as out:
        if len(RNAP_list) > 0:
            sigma = RNAP_list[-1].sigma
        else:
            sigma = (
                modelP.gene.Ldomain_lk - modelP.gene.Ldomain_lk_relaxed
            ) / modelP.gene.Ldomain_lk_relaxed
        out.write(
            "%s\t%.5f\t%d\n"
            % (
                str(traj.time),
                sigma,
                (len(RNAP_list) > 0 and RNAP_list[-1].t_elongating),
            )
        )
    return


def write_transcripts_on_the_fly(traj, simuP, header=False):
    """mean properties for every stage of the transcription process"""

    fi = simuP.fo_out + "/mean_properties.txt"

    if header:
        with open(fi, "w") as out:
            out.write(
                "time\ttranscripts_np\tprod_rate\tmean_prod_time\tmean_bind_time\tmean_ocf_time\tmean_esc_time\tmean_init_time\tmean_elong_time\n"
            )
        return

    with open(fi, "a") as out:
        prod_rate = traj.Ntranscripts / traj.time
        out.write(
            "%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
            % (
                traj.time,
                traj.Ntranscripts,
                prod_rate,
                traj.prod_times["mean"],
                traj.b2b_times["mean"],
                traj.ocf_times["mean"],
                traj.esc_times["mean"],
                traj.initiation_times["mean"],
                traj.elongation_times["mean"],
            )
        )

    return


def write_statistics(traj, simuP):
    """some statistics"""

    if traj.sigma_at_ocf:
        with open(simuP.fo_out + "/sigma_at_ocf_stage.txt", "w") as out:
            for s in traj.sigma_at_ocf:
                out.write("%f\n" % s)

    return


# I/O
def output_variables(cmd, modelP, simuP):
    """writing out parameters and variables"""
    
    with open(simuP.fo_out + "/param_var.txt", "w") as out:
        out.write(cmd + "\n")

        out.write("\n")
        out.write("###########\n")
        out.write("# General #\n")
        out.write("###########\n")
        dicto = modelP.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("########\n")
        out.write("# Gene #\n")
        out.write("########\n")
        dicto = modelP.gene.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("############\n")
        out.write("# Promoter #\n")
        out.write("############\n")
        dicto = modelP.promoter.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("########\n")
        out.write("# RNAP #\n")
        out.write("########\n")
        dicto = modelP.rnap.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("#########\n")
        out.write("# TopoI #\n")
        out.write("#########\n")
        dicto = modelP.topoI.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("##########\n")
        out.write("# Gyrase #\n")
        out.write("##########\n")
        dicto = modelP.gyrase.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("#######\n")
        out.write("# DNA #\n")
        out.write("#######\n")
        dicto = modelP.dna.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("###################\n")
        out.write("# Coarse graining #\n")
        out.write("###################\n")
        dicto = modelP.coarse_g.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")

        out.write("\n")
        out.write("##############\n")
        out.write("# Statistics #\n")
        out.write("##############\n")
        dicto = simuP.__dict__
        for key, val in dicto.items():
            out.write(key + "\t" + str(val) + "\n")
