import numpy as np, os, sys
import math
from param_var import *

"""
    Functions used by ../bin/twin.py
    For lightness, only downstream quantities are referred to as _down
    i.e., upstream quantities are NOT referred to as _up -- this might change in the future
"""

# supercoiling densities
def _sigma(rnap: RNAP, loc='up'):
    return (rnap.Lk[loc] - rnap.Lk0[loc]) / rnap.Lk0[loc]

# Simulation runs
def generate_run_follow_promoter(modelP: ModelParam, simuP: SimuParam):

    traj = Trajectory()
    RNAP_list = []
    # DNA-bound RNAPs

    nextevent2iter = {tag: -1 for tag in ["b", "oc", "esc"]}
    # next trial for binding, OC formation and promoter escape
    # -1 to avoid initial True
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = int(Z / modelP.coarse_g.tau_0)

    write_follow_promoter(traj, None, modelP, simuP, header=True)
    while traj.niter < simuP.Niterations and traj.Ntranscripts < simuP.Ntranscripts_max:
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


def generate_run_multiple_transcrtipts(modelP: ModelParam, simuP: SimuParam):

    traj = Trajectory()
    RNAP_list = []
    # DNA-bound RNAPs

    nextevent2iter = {tag: -1 for tag in ["b", "oc", "esc"]}
    # nextevent2iter: next trial for binding, OC formation and promoter escape
    # -1 to avoid initial True
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = int(Z / modelP.coarse_g.tau_0)

    write_transcripts_on_the_fly(traj, simuP, header=True)
    while traj.niter < simuP.Niterations and traj.Ntranscripts < simuP.Ntranscripts_max:
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

    return


# Transcription stages
def binding_stage(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter):
    """Binding of the RNAP at the promoter"""

    if (
        not RNAP_list
        or RNAP_list[-1].X
        > modelP.gene.rnap_xi + modelP.rnap.excluded_length_elongation
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
        rnap.Lk0['up'] = modelP.gene.Lk0_rnap_xi

        # SIGMA PROPERTIES
        if not RNAP_list:
            # a single RNAP => topological properties dictated by the entire topological domain
            rnap.sigma['up'] = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
            ) / modelP.gene.Lk0_domain
            rnap.Lk0['down'] = modelP.gene.Lk0_domain - modelP.gene.Lk0_rnap_xi
        else:
            # multiple RNAP => topological properties dictated by the immediate downstream RNAP
            rnap.sigma['up'] = RNAP_list[-1].sigma['up']
            rnap.Lk0['down'] = RNAP_list[-1].Lk0['up'] - rnap.Lk0['up']

        rnap.Lk['up'] = (1 + rnap.sigma['up']) * rnap.Lk0['up']

        rnap.sigma['down'] = rnap.sigma['up']
        rnap.Lk['down'] = (1 + rnap.sigma['down']) * rnap.Lk0['down']

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


def oc_formation_stage(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter):
    """OC formation if sigma <= threshold"""

    if RNAP_list[-1].sigma['up'] <= modelP.promoter.sigma_o:
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


def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
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
        RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
        RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']

    return


def topo_stage(modelP: ModelParam, RNAP_list):
    """TopoI and gyrase activity"""

    if RNAP_list:
        topo_stage_RNAPpresent(modelP, RNAP_list)
    else:
        topo_stage_RNAPabsent(modelP)

    return


def topo_stage_RNAPpresent(modelP: ModelParam, RNAP_list):
    """TopoI and gyrase activity in the presence of at least one DNA-bound RNAP"""

    # UPSTREAM
    # non-specific activities
    DTopoI, DGyrase = 0, 0
    if not RNAP_list[-1].t_elongating:
        # the most upstream RNAP (at the promoter) is not a barrier
        if len(RNAP_list) == 1:
            # topoisomerases can act anywhere along the domain
            domain_length_topo = modelP.gene.L_domain
            DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
            DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
        else:
            # the second RNAP is a barrier and we consider activity upstream
            domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
            DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
            # RNAP_list[-1] because RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up'] here
            DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
            # RNAP_list[-1] because RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up'] here
    else:
        # the most upstream RNAP is a barrier and we consider activity upstream
        domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
        DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
        DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)

    # specific activity in the presence of transcription (only TopoI)
    DTopoI_spec = 0
    if len(RNAP_list) > 1 or RNAP_list[-1].t_elongating:
        DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)

    if DTopoI != 0 or DGyrase != 0 or DTopoI_spec != 0:  # updating topo properties
        modelP.gene.Lk_domain += DTopoI + DGyrase + DTopoI_spec

        if not RNAP_list[-1].t_elongating:
            # properties of non-elongating RNAP are dictated by its downstream RNAP (if it exists)
            if len(RNAP_list) == 1:
                RNAP_list[0].sigma['up'] = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
            else:
                RNAP_list[-2].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
                RNAP_list[-2].sigma['up'] = _sigma(RNAP_list[-2], 'up')
                RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']

            RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']
            RNAP_list[-1].sigma['down'] = RNAP_list[
                -1
            ].sigma['up']  # because RNAP is not a barrier
            RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[
                -1
            ].Lk0['down']
        else:
            # RNAP is a barrier
            RNAP_list[-1].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
            RNAP_list[-1].sigma['up'] = _sigma(RNAP_list[-1], 'up')

    # DOWNSTREAM
    DTopoI_down, DGyrase_down, DGyrase_spec = 0, 0, 0
    if RNAP_list[0].t_elongating:
        # if non elongating, this means a single non-elongating RNAP => treated at the upstream level
        domain_length_topo = RNAP_list[0].Lk0['down'] * modelP.dna.n
        DTopoI_down = DLk_TopoI(domain_length_topo, RNAP_list[0], modelP, loc="down")
        DGyrase_down = DLk_Gyrase(domain_length_topo, RNAP_list[0], modelP, loc="down")
        DGyrase_spec = DLk_Gyrase("spec", RNAP_list[0], modelP, loc="down")

        modelP.gene.Lk_domain += DTopoI_down + DGyrase_down + DGyrase_spec
        RNAP_list[0].Lk['down'] += DTopoI_down + DGyrase_down + DGyrase_spec
        RNAP_list[0].sigma['down'] = _sigma(RNAP_list[0], 'down')

    return


def topo_stage_RNAPabsent(modelP: ModelParam):
    """TopoI and gyrase activity in the absence of RNAP"""

    # TopoI
    sigma = (
        modelP.gene.Lk_domain - modelP.gene.Lk0_domain
    ) / modelP.gene.Lk0_domain
    modelP.gene.Lk_domain += DLk_TopoI_noRNAP(modelP.gene.L_domain, modelP, sigma)

    # gyrase
    sigma = (
        modelP.gene.Lk_domain - modelP.gene.Lk0_domain
    ) / modelP.gene.Lk0_domain
    modelP.gene.Lk_domain += DLk_Gyrase_noRNAP(modelP.gene.L_domain, modelP, sigma)

    return


# Elementary generations of linking numbers
def DLk_TopoI(domain_length_topo, rnap: RNAP, modelP: ModelParam, loc="up"):
    """
    TopoI activity associated with an RNAP
    - worked for both upstream and downstream the "RNAP convoy"
    - not active if sigma > sigma_active
    """

    if rnap.sigma[loc] > modelP.topoI.sigma_active:
        return 0
    else:
        if domain_length_topo != "spec":
            return np.random.poisson(
                modelP.coarse_g.p_topoI_ns_per_bp * domain_length_topo
            )
        else:
            # idpt of distance
            return np.random.uniform() <= modelP.coarse_g.p_topoI_s


def DLk_Gyrase(domain_length_topo, rnap: RNAP, modelP: ModelParam, loc="up"):
    """
    Gyrase activity associated with an RNAP
    - worked for both upstream and downstream the "RNAP convoy"
    - not active if sigma < sigma_stall
    """

    if rnap.sigma[loc] < modelP.rnap.sigma_stall:
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


def DLk_TopoI_noRNAP(domain_length_topo, modelP: ModelParam, sigma):
    """
    TopoI activity in the absence of any RNAP => sigma is specified as an argument
    - not active if sigma > sigma_active
    """

    if sigma > modelP.topoI.sigma_active:
        return 0
    else:
        return np.random.poisson(modelP.coarse_g.p_topoI_ns_per_bp * domain_length_topo)


def DLk_Gyrase_noRNAP(domain_length_topo, modelP: ModelParam, sigma):
    """
    Gyrase activity in the absence of any RNAP => sigma is specified as an argument
    - not active if sigma > sigma_active
    """

    if sigma < modelP.rnap.sigma_stall:
        # stalling torque is the gyrase threshold (not 0, otherwsie, we may be stuck if we start from sigma = 0)
        return 0
    else:
        return -2 * np.random.poisson(
            modelP.coarse_g.p_gyrase_ns_per_bp * domain_length_topo
        )


def elongation_stage(modelP: ModelParam, RNAP_list):
    """RNAPs translocation"""

    if RNAP_list:
        ix_ = np.arange(len(RNAP_list))
        np.random.shuffle(ix_)

        for ix_rnap in ix_:
            RNAP_translocation(ix_rnap, RNAP_list, modelP)

    return


# Translocations and their topological consequences
def RNAP_translocation(ix_rnap, RNAP_list, modelP: ModelParam):

    if (
        not RNAP_list[ix_rnap].t_elongating
        or not RNAP_list[ix_rnap].sigma['up'] >= modelP.rnap.sigma_stall
        or not RNAP_list[ix_rnap].sigma['down'] <= np.abs(modelP.rnap.sigma_stall)
    ):
        # not elongating or beyond sigma_stall = no translocation
        return

    # Rem1: linking numbers do not change, supercoiling density does
    # Rem2: no constrain on one RNAP overtaking another one (supercoiling constraints do the job)

    RNAP_list[ix_rnap].X += modelP.coarse_g.dx
    RNAP_list[ix_rnap].Lk0['up'] += modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma['up'] = _sigma(RNAP_list[ix_rnap], 'up')

    RNAP_list[ix_rnap].Lk0['down'] -= modelP.coarse_g.dLk
    RNAP_list[ix_rnap].sigma['down'] = _sigma(RNAP_list[ix_rnap], 'down')

    # UPDATING DOWNSTREAM RNAP
    if ix_rnap > 0:
        RNAP_list[ix_rnap - 1].Lk0['up'] -= modelP.coarse_g.dLk
        RNAP_list[ix_rnap - 1].sigma['up'] = _sigma(RNAP_list[ix_rnap - 1], 'up')

    # UPDATING UPSTREAM RNAP
    if ix_rnap < len(RNAP_list) - 1:
        if RNAP_list[ix_rnap + 1].t_elongating:
            # the RNAP is elongating
            RNAP_list[ix_rnap + 1].Lk0['down'] += modelP.coarse_g.dLk
            RNAP_list[ix_rnap + 1].sigma['down'] = _sigma(
                RNAP_list[ix_rnap + 1], 'down'
            )
        elif (
            ix_rnap == len(RNAP_list) - 2
            and not RNAP_list[len(RNAP_list) - 1].t_elongating
        ):
            # the RNAP is NON-elongating
            RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']
            RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']

            RNAP_list[-1].Lk0['down'] += modelP.coarse_g.dLk
            RNAP_list[-1].sigma['down'] = RNAP_list[-2].sigma['up']
            RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[
                -1
            ].Lk0['down']

    return


def termination_stage(modelP: ModelParam, simuP: SimuParam, RNAP_list, traj: Trajectory):
    """Termination stage: transcript production by the most downstream RNAP"""

    if RNAP_list and RNAP_list[0].X >= modelP.gene.term:
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

        # TRANSCRIPT TERMINATION
        # 1. We update the upstream RNAP, if it exists        
        if len(RNAP_list) > 1: # at least 2 RNAPs: upating of the upstream RNAP (index = 1)
            RNAP_list[1].Lk0['down'] += RNAP_list[0].Lk0['down']
            if RNAP_list[1].t_elongating:
                # only downstream properties are updated                
                RNAP_list[1].Lk['down'] += RNAP_list[0].Lk['down']
                RNAP_list[1].sigma['down'] = _sigma(RNAP_list[1], 'down')
            else:
                # both upstrean and dowsntream properties are updtaed from sigma of the domain
                RNAP_list[1].sigma['up'] = (
                    modelP.gene.Lk_domain - modelP.gene.Lk0_domain
                ) / modelP.gene.Lk0_domain
                RNAP_list[1].Lk['up'] = (1 + RNAP_list[1].sigma['up']) * RNAP_list[1].Lk0['up']

                RNAP_list[1].sigma['down'] = RNAP_list[1].sigma['up']
                RNAP_list[1].Lk['down'] = (1 + RNAP_list[1].sigma['down']) * RNAP_list[1].Lk0['down']
        # 2: We remove the RNAP
        del RNAP_list[0]

        if not traj.Ntranscripts % simuP.Nevery_transcripts:
            write_transcripts_on_the_fly(traj, simuP)

    return


# I/O
def verbosing(simuP: SimuParam, traj: Trajectory, RNAP_list):
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
                    RNAP_list[-1].sigma['up'],
                    RNAP_list[-1].sigma['down'],
                    end="\r",
                )
            else:
                print(
                    traj.niter,
                    traj.Ntranscripts,
                    len(RNAP_list),
                    RNAP_list[-1].sigma['up'],
                    RNAP_list[-1].sigma['down'],
                    RNAP_list[-1].X,
                    RNAP_list[-1].t_elongating,
                    RNAP_list[-2].sigma['up'],
                    RNAP_list[-2].sigma['down'],
                    RNAP_list[-2].X,
                    RNAP_list[-2].t_elongating,
                    RNAP_list[0].sigma['down'],
                    RNAP_list[0].X,
                    end="\r",
                )


def write_follow_promoter(traj: Trajectory, RNAP_list, modelP: ModelParam, simuP: SimuParam, header=False):
    """promoter properties"""

    fi = simuP.fo_out + "/traj_promoter.txt"

    if header:
        with open(fi, "w") as out:
            out.write("time\tsigma\tt_upRNAPelong\n")
        return

    with open(fi, "a") as out:
        if len(RNAP_list) > 0:
            sigma = RNAP_list[-1].sigma['up']
        else:
            sigma = (
                modelP.gene.Lk_domain - modelP.gene.Lk0_domain
            ) / modelP.gene.Lk0_domain
        out.write(
            "%s\t%.5f\t%d\n"
            % (
                str(traj.time),
                sigma,
                (len(RNAP_list) > 0 and RNAP_list[-1].t_elongating),
            )
        )
    return


def write_transcripts_on_the_fly(traj: Trajectory, simuP: SimuParam, header=False):
    """mean properties for every stage of the transcription process"""

    fi = simuP.fo_out + "/mean_properties.txt"

    if header:
        with open(fi, "w") as out:
            out.write(
                "transcripts_nb\ttime\tprod_rate\tmean_prod_time\tmean_bind_time\tmean_ocf_time\tmean_esc_time\tmean_init_time\tmean_elong_time\n"
            )
        return

    with open(fi, "a") as out:
        prod_rate = traj.Ntranscripts / traj.time
        out.write(
            "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"
            % (
                traj.Ntranscripts,
                traj.time,
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

# I/O
def output_variables(cmd, modelP: ModelParam, simuP: SimuParam):
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
