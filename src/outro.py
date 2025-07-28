import numpy as np, os, sys
import math
from param_var import *

"""
    Functions used by ../bin/twin.py    
"""
_stats = {
    'bind_oc': 0,
    'bind_bind': 0,
    'spec_events': 0,
    'unbind_events': 0,
    'iter': 0
}

def _dump_stats(label):
    print(f"\n=== STATS {label} ===")
    for k,v in _stats.items():
        print(f"{k}: {v}")

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
        
        # NOVO: TopoI unbinding
        topo1_unbinding_stage(modelP, RNAP_list)

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

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # Binding
        if traj.niter == nextevent2iter["b"]:
            binding_stage(modelP, RNAP_list, traj, nextevent2iter)

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # OC formation
        if traj.niter == nextevent2iter["oc"]:
            oc_formation_stage(modelP, RNAP_list, traj, nextevent2iter)

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # Promoter escape
        if traj.niter == nextevent2iter["esc"]:
            escape_stage(modelP, RNAP_list, traj)

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # Topoisomerases
        topo_stage(modelP, RNAP_list)

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # Elongation
        elongation_stage(modelP, RNAP_list)

        # NOVO: TopoI unbinding
        topo1_unbinding_stage(modelP, RNAP_list)

        # Termination
        termination_stage(modelP, simuP, RNAP_list, traj)

        # NOVO: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        traj.niter += 1
        traj.time = traj.niter * modelP.coarse_g.tau_0

    # NOVO: TopoI unbinding
    #topo1_unbinding_stage(modelP, RNAP_list)

    _dump_stats("multiple_transcripts")
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
        #rnap.topo1_bound = True #novo
        #rnap.topo_spec.bind()
        _stats['bind_bind'] += 1

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
        rnap = RNAP_list[-1]
        #rnap.topo1_bound = True
        rnap.topo_spec.bind()
        _stats['bind_oc'] += 1
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


# def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
#     """promoter escape => RNAP is now in elongating mode"""

#     RNAP_list[-1].t_elongating = True
#     RNAP_list[-1].tesc = traj.time

#     traj.esc_times["mean"] = (
#         traj.esc_times["n"] * traj.esc_times["mean"]
#         + RNAP_list[-1].tesc
#         - RNAP_list[-1].tocf
#     ) / (traj.esc_times["n"] + 1)
#     traj.esc_times["n"] += 1

#     traj.initiation_times["mean"] = (
#         traj.initiation_times["n"] * traj.initiation_times["mean"]
#         + RNAP_list[-1].tesc
#         - RNAP_list[-1].tb
#     ) / (traj.initiation_times["n"] + 1)
#     traj.initiation_times["n"] += 1

#     # UPDATING DOWNSTREAM RNAP
#     # change Lk because the new elongating RNA becomes a barrier
#     # one can actually check the conservation of the Lk
#     if len(RNAP_list) > 1:
#         RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
#         RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']

#     return


def topo_stage(modelP: ModelParam, RNAP_list):
    """TopoI and gyrase activity"""

    if RNAP_list:
        topo_stage_RNAPpresent(modelP, RNAP_list)
    else:
        topo_stage_RNAPabsent(modelP)

    return


# def topo_stage_RNAPpresent(modelP: ModelParam, RNAP_list):
#     """TopoI and gyrase activity in the presence of at least one DNA-bound RNAP"""

#     # UPSTREAM
#     # non-specific activities
#     DTopoI, DGyrase = 0, 0
#     if not RNAP_list[-1].t_elongating:
#         # the most upstream RNAP (at the promoter) is not a barrier
#         if len(RNAP_list) == 1:
#             # topoisomerases can act anywhere along the domain
#             domain_length_topo = modelP.gene.L_domain
#             DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#             DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
#         else:
#             # the second RNAP is a barrier and we consider activity upstream
#             domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#             DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#             # RNAP_list[-1] because RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up'] here
#             DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
#             # RNAP_list[-1] because RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up'] here
#     else:
#         # the most upstream RNAP is a barrier and we consider activity upstream
#         domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
#         DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)

#     # specific activity in the presence of transcription (only TopoI)
#     DTopoI_spec = 0
#     if len(RNAP_list) > 1 or RNAP_list[-1].t_elongating:
#         DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)

#     if DTopoI != 0 or DGyrase != 0 or DTopoI_spec != 0:  # updating topo properties
#         modelP.gene.Lk_domain += DTopoI + DGyrase + DTopoI_spec

#         if not RNAP_list[-1].t_elongating:
#             # properties of non-elongating RNAP are dictated by its downstream RNAP (if it exists)
#             if len(RNAP_list) == 1:
#                 RNAP_list[0].sigma['up'] = (
#                     modelP.gene.Lk_domain - modelP.gene.Lk0_domain
#                 ) / modelP.gene.Lk0_domain
#             else:
#                 RNAP_list[-2].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
#                 RNAP_list[-2].sigma['up'] = _sigma(RNAP_list[-2], 'up')
#                 RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']

#             RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']
#             RNAP_list[-1].sigma['down'] = RNAP_list[
#                 -1
#             ].sigma['up']  # because RNAP is not a barrier
#             RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[
#                 -1
#             ].Lk0['down']
#         else:
#             # RNAP is a barrier
#             RNAP_list[-1].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
#             RNAP_list[-1].sigma['up'] = _sigma(RNAP_list[-1], 'up')

#     # DOWNSTREAM
#     DTopoI_down, DGyrase_down, DGyrase_spec = 0, 0, 0
#     if RNAP_list[0].t_elongating:
#         # if non elongating, this means a single non-elongating RNAP => treated at the upstream level
#         domain_length_topo = RNAP_list[0].Lk0['down'] * modelP.dna.n
#         DTopoI_down = DLk_TopoI(domain_length_topo, RNAP_list[0], modelP, loc="down")
#         DGyrase_down = DLk_Gyrase(domain_length_topo, RNAP_list[0], modelP, loc="down")
#         DGyrase_spec = DLk_Gyrase("spec", RNAP_list[0], modelP, loc="down")

#         modelP.gene.Lk_domain += DTopoI_down + DGyrase_down + DGyrase_spec
#         RNAP_list[0].Lk['down'] += DTopoI_down + DGyrase_down + DGyrase_spec
#         RNAP_list[0].sigma['down'] = _sigma(RNAP_list[0], 'down')

#     return


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






# def update_presence_topo1_spec(rnap, modelP, dt):
#     """
#     Atualiza estocasticamente rnap.topo1_bound:
#       - rnap.topo1_bound == 1 enquanto a TopoI permanecer ligada
#       - com probabilidade k_unbind * dt ela se desliga (volta para 0)
#     """
#     # garante que o atributo existe
#     if not hasattr(rnap, 'topo1_bound'):
#         rnap.topo1_bound = 1  # se não existia, assumimos ligada inicialmente

#     if rnap.topo1_bound == 1:
#         p_unbind = modelP.topoI.k_unbind * dt
#         # rola o dado
#         if np.random.rand() < p_unbind:
#             rnap.topo1_bound = 0
#             # opcional: log para debug
#             # print(f"[UNBIND] t={modelP.current_time:.2f}  topo1_bound→0")


# def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
#     """promoter escape => RNAP entra em elongação e TopoI específico se liga"""
#     # original :contentReference[oaicite:0]{index=0}
#     RNAP_list[-1].t_elongating = True
#     RNAP_list[-1].tesc = traj.time

#     # estatísticas de escape e initiation
#     traj.esc_times["mean"] = (
#         traj.esc_times["n"] * traj.esc_times["mean"]
#         + RNAP_list[-1].tesc - RNAP_list[-1].tocf
#     ) / (traj.esc_times["n"] + 1)
#     traj.esc_times["n"] += 1

#     traj.initiation_times["mean"] = (
#         traj.initiation_times["n"] * traj.initiation_times["mean"]
#         + RNAP_list[-1].tesc - RNAP_list[-1].tb
#     ) / (traj.initiation_times["n"] + 1)
#     traj.initiation_times["n"] += 1

#     # agora que passou ao modo elongação, ligamos TopoI específico
#     RNAP_list[-1].topo1_bound = 1

#     # UPDATING DOWNSTREAM RNAP (mantém o Lk)
#     if len(RNAP_list) > 1:
#         RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
#         RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']


# def topo_stage_RNAPpresent(modelP: ModelParam, RNAP_list):
#     """
#     TopoI e gyrase em presença de RNAP(s) ligados ao DNA,
#     agora com k_unbinding de TopoI específico.
#     """
#     # introduz update estocástico de unbinding
#     dt = modelP.coarse_g.tau_0
#     for rnap in RNAP_list:
#         update_presence_topo1_spec(rnap, modelP, dt)

#     # UPSTREAM não-específico
#     DTopoI, DGyrase = 0, 0
#     if not RNAP_list[-1].t_elongating:
#         if len(RNAP_list) == 1:
#             domain_length = modelP.gene.L_domain
#         else:
#             domain_length = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#     else:
#         domain_length = RNAP_list[-1].Lk0['up'] * modelP.dna.n

#     DTopoI  = DLk_TopoI  (domain_length, RNAP_list[-1], modelP)
#     DGyrase = DLk_Gyrase (domain_length, RNAP_list[-1], modelP)

#     # específico (só TopoI), só se TopoI ainda estiver ligada
#     DTopoI_spec = 0
#     last = RNAP_list[-1]
#     if (len(RNAP_list) > 1 or last.t_elongating) and getattr(last, 'topo1_bound', 0) == 1:
#         DTopoI_spec = DLk_TopoI("spec", last, modelP)  # branch “spec” em :contentReference[oaicite:1]{index=1}

#     # aplica todas as mudanças de linking
#     if (DTopoI + DGyrase + DTopoI_spec) != 0:
#         modelP.gene.Lk_domain += DTopoI + DGyrase + DTopoI_spec

#         if not last.t_elongating:
#             if len(RNAP_list) == 1:
#                 last.sigma['up'] = (modelP.gene.Lk_domain - modelP.gene.Lk0_domain) / modelP.gene.Lk0_domain
#             else:
#                 prev = RNAP_list[-2]
#                 prev.Lk['up'] += DTopoI + DGyrase + DTopoI_spec
#                 prev.sigma['up'] = _sigma(prev, 'up')
#                 last.sigma['up'] = prev.sigma['up']

#             last.Lk['up']   = (1 + last.sigma['up'])   * last.Lk0['up']
#             last.sigma['down'] = last.sigma['up']
#             last.Lk['down'] = (1 + last.sigma['down']) * last.Lk0['down']
#         else:
#             last.Lk['up']   += DTopoI + DGyrase + DTopoI_spec
#             last.sigma['up'] = _sigma(last, 'up')

#     # DOWNSTREAM (idem original)
#     D_down, G_down, G_sp = 0, 0, 0
#     first = RNAP_list[0]
#     if first.t_elongating:
#         length_down = first.Lk0['down'] * modelP.dna.n
#         D_down  = DLk_TopoI  (length_down, first, modelP, loc="down")
#         G_down  = DLk_Gyrase (length_down, first, modelP, loc="down")
#         G_sp    = DLk_Gyrase ("spec", first, modelP, loc="down")

#         modelP.gene.Lk_domain += D_down + G_down + G_sp
#         first.Lk['down']      += D_down + G_down + G_sp
#         first.sigma['down']    = _sigma(first, 'down')



# def update_presence_topo1_spec(rnap, modelP, dt): NOVO NOVO
#     """
#     Atualiza estocasticamente rnap.topo1_bound:
#       - rnap.topo1_bound == 1 enquanto a TopoI permanecer ligada
#       - com probabilidade k_unbind * dt ela se desliga (ponte para 0)
#     """
#     # garante que o atributo existe
#     if not hasattr(rnap, "topo1_bound"):
#         rnap.topo1_bound = 0

#     # se estiver ligada, tenta descolar
#     if rnap.topo1_bound == 1:
#         if np.random.rand() < modelP.topoI.k_unbind * dt:
#             rnap.topo1_bound = 0



# def escape_stage(modelP, RNAP_list, traj):
#     """
#     Promoter escape:
#       - o último RNAP da lista sai do promoter e vira elongating
#       - registra tempo de escape em traj.esc_times
#       - dispara binding imediato de TopoI específico
#     """

#     # ————————————————————————————————————————————————
#     # 0) Inicializa o dicionário de estatísticas de escape, se preciso
#     if not hasattr(traj, "esc_times") or not isinstance(traj.esc_times, dict):
#         traj.esc_times = {"sum": 0.0, "count": 0, "mean": 0.0}
#     else:
#         # garante as chaves mínimas
#         for k in ("sum", "count", "mean"):
#             traj.esc_times.setdefault(k, 0.0 if k == "mean" or k == "sum" else 0)

#     # ————————————————————————————————————————————————
#     # 1) Promoter escape: o último RNAP começa a elongar
#     rnap = RNAP_list[-1]
#     rnap.t_elongating = True
#     rnap.tesc       = traj.time

#     # 2) Atualiza estatísticas de tempo de escape
#     traj.esc_times["sum"]   += rnap.tesc
#     traj.esc_times["count"] += 1
#     traj.esc_times["mean"]   = traj.esc_times["sum"] / traj.esc_times["count"]

#     # ————————————————————————————————————————————————
#     # 3) Novo: TopoI específico se liga imediatamente a este RNAP
#     #     (será liberada estocasticamente em topo_stage_RNAPpresent)
#     rnap.topo1_bound = 1

#     # nada é retornado; o estado de traj e RNAP_list é modificado “in place”


# def topo_stage_RNAPpresent(modelP, RNAP_list):
#     """
#     TopoI e Gyrase atuam enquanto há pelo menos um RNAP ligado ao DNA.
#     Inclui agora unbinding estocástico de TopoI específico.
#     """
#     # ─── 0) DESEMBARALHA TopoI específico (unbinding) ───
#     dt = modelP.coarse_g.tau_0
#     for rnap in RNAP_list:
#         update_presence_topo1_spec(rnap, modelP, dt)

#     # ─── I) UPSTREAM — atividades não‐específicas ───
#     DTopoI, DGyrase = 0, 0
#     last = RNAP_list[-1]
#     if not last.t_elongating:
#         if len(RNAP_list) == 1:
#             domain_length_topo = modelP.gene.L_domain
#         else:
#             domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#         DTopoI  = DLk_TopoI(domain_length_topo, last, modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, last, modelP)
#     else:
#         domain_length_topo = last.Lk0['up'] * modelP.dna.n
#         DTopoI  = DLk_TopoI(domain_length_topo, last, modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, last, modelP)

#     # ─── II) UPSTREAM — atividade específica de TopoI ───
#     DTopoI_spec = 0
#     # só age se houver algum RNAP em elongação _e_ TopoI ainda ligada
#     if (len(RNAP_list) > 1 or last.t_elongating) and getattr(last, "topo1_bound", 0) == 1:
#         DTopoI_spec = DLk_TopoI("spec", last, modelP)

#     # ─── III) AGORA, atualiza torção global e local ───
#     if DTopoI or DGyrase or DTopoI_spec:
#         modelP.gene.Lk_domain += DTopoI + DGyrase + DTopoI_spec

#         # propaga para os RNAPs conforme elongação ou não
#         if not last.t_elongating:
#             if len(RNAP_list) == 1:
#                 # único RNAP não‐elongando
#                 RNAP_list[0].sigma['up'] = (
#                     modelP.gene.Lk_domain - modelP.gene.Lk0_domain
#                 ) / modelP.gene.Lk0_domain
#             else:
#                 prev = RNAP_list[-2]
#                 prev.Lk['up']   += DTopoI + DGyrase + DTopoI_spec
#                 prev.sigma['up'] = _sigma(prev, 'up')
#                 last.sigma['up'] = prev.sigma['up']

#             last.Lk['up']    = (1 + last.sigma['up']) * last.Lk0['up']
#             last.sigma['down'] = last.sigma['up']
#             last.Lk['down']  = (1 + last.sigma['down']) * last.Lk0['down']
#         else:
#             # barrier de torção
#             last.Lk['up']   += DTopoI + DGyrase + DTopoI_spec
#             last.sigma['up'] = _sigma(last, 'up')

#     # ─── IV) DOWNSTREAM — apenas se o RNAP estiver elongando ───
#     first = RNAP_list[0]
#     if first.t_elongating:
#         domain_length_topo = first.Lk0['down'] * modelP.dna.n
#         D_down  = DLk_TopoI(domain_length_topo, first, modelP, loc="down")
#         G_down  = DLk_Gyrase(domain_length_topo, first, modelP, loc="down")
#         G_sp    = DLk_Gyrase("spec", first, modelP, loc="down")

#         modelP.gene.Lk_domain += D_down + G_down + G_sp
#         first.Lk['down']    += D_down + G_down + G_sp
#         first.sigma['down']  = _sigma(first, 'down')

#     # não retorna nada


# def _ensure_topo1_bound(rnap):  MODELO BOM ESTA QUASE
#     """Função auxiliar para garantir que o atributo topo1_bound existe"""
#     if not hasattr(rnap, "topo1_bound"):
#         rnap.topo1_bound = 0

# def update_presence_topo1_spec(rnap, modelP, dt):
#     """
#     Atualiza estocasticamente rnap.topo1_bound:
#       - rnap.topo1_bound == 1 enquanto a TopoI permanecer ligada
#       - com probabilidade k_unbind * dt ela se desliga (ponte para 0)
#     """
#     # Garante que o atributo existe
#     _ensure_topo1_bound(rnap)

#     # Se estiver ligada, tenta desligar
#     if rnap.topo1_bound == 1:
#         # Verifica se dt não é muito pequeno para evitar problemas numéricos
#         unbind_prob = modelP.topoI.k_unbind * dt
#         if unbind_prob > 1.0:
#             # Se a probabilidade for > 1, força o unbinding
#             rnap.topo1_bound = 0
#         elif np.random.rand() < unbind_prob:
#             rnap.topo1_bound = 0

# def escape_stage(modelP, RNAP_list, traj):
#     """
#     Promoter escape:
#       - o último RNAP da lista sai do promoter e vira elongating
#       - registra tempo de escape em traj.esc_times
#       - dispara binding imediato de TopoI específico
#     """
    
#     # ————————————————————————————————————————————————
#     # 0) Verifica se há RNAPs na lista
#     if not RNAP_list:
#         return
    
#     # ————————————————————————————————————————————————
#     # 1) Inicializa as estatísticas de escape se necessário
#     if not hasattr(traj, "esc_times"):
#         traj.esc_times = {"sum": 0.0, "count": 0, "mean": 0.0}
    
#     # Versão compatível com o código original que usa "n" em vez de "count"
#     if "n" in traj.esc_times:
#         # Mantém compatibilidade com estrutura original
#         traj.esc_times["mean"] = (
#             traj.esc_times["n"] * traj.esc_times["mean"]
#             + traj.time - RNAP_list[-1].tocf
#         ) / (traj.esc_times["n"] + 1)
#         traj.esc_times["n"] += 1
#     else:
#         # Nova estrutura
#         for k in ("sum", "count", "mean"):
#             traj.esc_times.setdefault(k, 0.0 if k in ("mean", "sum") else 0)

#     # ————————————————————————————————————————————————
#     # 2) Promoter escape: o último RNAP começa a elongar
#     rnap = RNAP_list[-1]
#     rnap.t_elongating = True
#     rnap.tesc = traj.time

#     # 3) Atualiza estatísticas de tempo de escape (nova estrutura)
#     if "n" not in traj.esc_times:
#         traj.esc_times["sum"] += rnap.tesc - rnap.tocf
#         traj.esc_times["count"] += 1
#         traj.esc_times["mean"] = traj.esc_times["sum"] / traj.esc_times["count"]

#     # ————————————————————————————————————————————————
#     # 4) Atualiza estatísticas de tempo de iniciação (compatibilidade)
#     if hasattr(traj, "initiation_times") and "n" in traj.initiation_times:
#         traj.initiation_times["mean"] = (
#             traj.initiation_times["n"] * traj.initiation_times["mean"]
#             + rnap.tesc - rnap.tb
#         ) / (traj.initiation_times["n"] + 1)
#         traj.initiation_times["n"] += 1

#     # ————————————————————————————————————————————————
#     # 5) TopoI específico se liga imediatamente a este RNAP
#     rnap.topo1_bound = 1

#     # ————————————————————————————————————————————————
#     # 6) Atualiza RNAP downstream (código original)
#     if len(RNAP_list) > 1:
#         RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
#         RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']

# def topo_stage_RNAPpresent(modelP, RNAP_list):
#     """
#     TopoI e Gyrase atuam enquanto há pelo menos um RNAP ligado ao DNA.
#     Inclui agora unbinding estocástico de TopoI específico.
#     """
    
#     # Verifica se há RNAPs
#     if not RNAP_list:
#         return
    
#     # ─── 0) UNBINDING de TopoI específico ───
#     dt = modelP.coarse_g.tau_0
#     for rnap in RNAP_list:
#         update_presence_topo1_spec(rnap, modelP, dt)

#     # ─── I) UPSTREAM — atividades não-específicas ───
#     DTopoI, DGyrase = 0, 0
#     last = RNAP_list[-1]
    
#     if not last.t_elongating:
#         if len(RNAP_list) == 1:
#             domain_length_topo = modelP.gene.L_domain
#         else:
#             domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#         DTopoI = DLk_TopoI(domain_length_topo, last, modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, last, modelP)
#     else:
#         domain_length_topo = last.Lk0['up'] * modelP.dna.n
#         DTopoI = DLk_TopoI(domain_length_topo, last, modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, last, modelP)

#     # ─── II) UPSTREAM — atividade específica de TopoI ───
#     DTopoI_spec = 0
#     # Só age se houver transcriçao ativa E TopoI ainda ligada
#     if (len(RNAP_list) > 1 or last.t_elongating):
#         # Garante que o atributo existe antes de verificar
#         _ensure_topo1_bound(last)
#         if last.topo1_bound == 1:
#             DTopoI_spec = DLk_TopoI("spec", last, modelP)

#     # ─── III) Atualiza torção global e local ───
#     total_change = DTopoI + DGyrase + DTopoI_spec
#     if total_change != 0:
#         modelP.gene.Lk_domain += total_change

#         # Propaga para os RNAPs conforme elongação ou não
#         if not last.t_elongating:
#             if len(RNAP_list) == 1:
#                 # Único RNAP não-elongando
#                 RNAP_list[0].sigma['up'] = (
#                     modelP.gene.Lk_domain - modelP.gene.Lk0_domain
#                 ) / modelP.gene.Lk0_domain
#             else:
#                 prev = RNAP_list[-2]
#                 prev.Lk['up'] += total_change
#                 prev.sigma['up'] = _sigma(prev, 'up')
#                 last.sigma['up'] = prev.sigma['up']

#             last.Lk['up'] = (1 + last.sigma['up']) * last.Lk0['up']
#             last.sigma['down'] = last.sigma['up']
#             last.Lk['down'] = (1 + last.sigma['down']) * last.Lk0['down']
#         else:
#             # Barreira de torção
#             last.Lk['up'] += total_change
#             last.sigma['up'] = _sigma(last, 'up')

#     # ─── IV) DOWNSTREAM — apenas se o RNAP estiver elongando ───
#     first = RNAP_list[0]
#     if first.t_elongating:
#         domain_length_topo = first.Lk0['down'] * modelP.dna.n
#         D_down = DLk_TopoI(domain_length_topo, first, modelP, loc="down")
#         G_down = DLk_Gyrase(domain_length_topo, first, modelP, loc="down")
#         G_sp = DLk_Gyrase("spec", first, modelP, loc="down")

#         downstream_change = D_down + G_down + G_sp
#         modelP.gene.Lk_domain += downstream_change
#         first.Lk['down'] += downstream_change
#         first.sigma['down'] = _sigma(first, 'down')


# def update_presence_topo1_spec(rnap, modelP, dt):
#     """
#     Atualiza estocasticamente rnap.topo1_bound:
#       - rnap.topo1_bound == 1 enquanto a TopoI permanecer ligada
#       - com probabilidade k_unbind * dt ela se desliga (→ 0)
#     Se k_unbind == 0 então unbind_prob == 0 e nunca desliga, reproduzindo o comportamento antigo.
#     """
#     # inicializa atributo, se precisar
#     if not hasattr(rnap, "topo1_bound"):
#         rnap.topo1_bound = 0

#     # se estiver ligada, testamos unbinding
#     if rnap.topo1_bound == 1:
#         p_unbind = modelP.topoI.k_unbind * dt
#         # probabilidade maior que 1 => força unbind
#         if p_unbind >= 1.0 or np.random.rand() < p_unbind:
#             rnap.topo1_bound = 0


# def escape_stage(modelP, RNAP_list, traj):
#     """
#     Promoter escape:
#       - o último RNAP sai do promoter e vira elongating
#       - registra tempo de escape em traj.esc_times (sum/count/mean)
#       - dispara binding imediato de TopoI específico
#       - ajusta Lk0 upstream, se houver
#     """
#     if not RNAP_list:
#         return

#     # ─── 1) Normaliza traj.esc_times para o formato {sum, count, mean} ───
#     if not hasattr(traj, "esc_times"):
#         # nunca inicializado: cria do zero
#         traj.esc_times = {"sum": 0.0, "count": 0, "mean": 0.0}
#     else:
#         et = traj.esc_times
#         if "n" in et:
#             # migrar do formato antigo (n/mean) → (sum/count/mean)
#             n = et["n"]
#             mean_old = et.get("mean", 0.0)
#             traj.esc_times = {
#                 "sum": n * mean_old,
#                 "count": n,
#                 "mean": mean_old,
#             }
#         else:
#             # garante chaves sum/count/mean
#             traj.esc_times.setdefault("sum", 0.0)
#             traj.esc_times.setdefault("count", 0)
#             traj.esc_times.setdefault("mean", 0.0)

#     # ─── 2) Computa intervalo de escape ───
#     rnap = RNAP_list[-1]
#     dt_escape = traj.time - rnap.tocf

#     # ─── 3) Atualiza sum/count/mean ───
#     traj.esc_times["sum"] += dt_escape
#     traj.esc_times["count"] += 1
#     traj.esc_times["mean"] = traj.esc_times["sum"] / traj.esc_times["count"]

#     # ─── 4) Marca como elongating e registra tesc ───
#     rnap.t_elongating = True
#     rnap.tesc = traj.time

#     # ─── 5) Atualiza estatísticas de iniciação (se existirem) ───
#     if hasattr(traj, "initiation_times"):
#         it = traj.initiation_times
#         # migra caso venha do formato antigo
#         if "n" in it:
#             n0 = it["n"]
#             m0 = it.get("mean", 0.0)
#             it.clear()
#             it["sum"] = n0 * m0
#             it["count"] = n0
#             it["mean"] = m0
#         it.setdefault("sum", 0.0)
#         it.setdefault("count", 0)
#         it.setdefault("mean", 0.0)

#         dt_init = rnap.tesc - rnap.tb
#         it["sum"] += dt_init
#         it["count"] += 1
#         it["mean"] = it["sum"] / it["count"]

#     # ─── 6) Binding imediato de TopoI específico ───
#     rnap.topo1_bound = 1

#     # ─── 7) Ajusta Lk0 upstream (se houver um RNAP anterior) ───
#     if len(RNAP_list) > 1:
#         prev = RNAP_list[-2]
#         prev.Lk0['up'] -= rnap.Lk0['up']
#         prev.Lk['up'] = (1 + prev.sigma['up']) * prev.Lk0['up']



# def topo_stage_RNAPpresent(modelP, RNAP_list):
#     """
#     Enquanto houver RNAPs ligados ao DNA, TopoI e Gyrase atuam:
#       - primeiro roda unbinding estocástico de TopoI específico
#       - depois calcula atividade não-específica (Upstream e Downstream)
#       - finalmente atividade específica de TopoI (se topo1_bound == 1)
#     """
#     if not RNAP_list:
#         return

#     # 0) unbinding de qualquer topo1 já ligado
#     dt = modelP.coarse_g.tau_0
#     for rnap in RNAP_list:
#         update_presence_topo1_spec(rnap, modelP, dt)

#     last = RNAP_list[-1]

#     # I) atividade não-específica upstream
#     if not last.t_elongating:
#         if len(RNAP_list) == 1:
#             L_up = modelP.gene.L_domain
#         else:
#             L_up = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#     else:
#         L_up = last.Lk0['up'] * modelP.dna.n

#     D_ns = DLk_TopoI(L_up, last, modelP)
#     G_ns = DLk_Gyrase(L_up, last, modelP)

#     # II) atividade específica (apenas TopoI) upstream
#     D_sp = 0
#     if (len(RNAP_list) > 1 or last.t_elongating) and getattr(last, "topo1_bound", 0) == 1:
#         D_sp = DLk_TopoI("spec", last, modelP)

#     # III) aplica mudanças upstream
#     delta_up = D_ns + G_ns + D_sp
#     if delta_up != 0:
#         modelP.gene.Lk_domain += delta_up

#         if not last.t_elongating:
#             if len(RNAP_list) == 1:
#                 last.sigma['up'] = (modelP.gene.Lk_domain - modelP.gene.Lk0_domain) / modelP.gene.Lk0_domain
#             else:
#                 prev = RNAP_list[-2]
#                 prev.Lk['up'] += delta_up
#                 prev.sigma['up'] = _sigma(prev, 'up')
#                 last.sigma['up'] = prev.sigma['up']

#             last.Lk['up'] = (1 + last.sigma['up']) * last.Lk0['up']
#             last.sigma['down'] = last.sigma['up']
#             last.Lk['down'] = (1 + last.sigma['down']) * last.Lk0['down']
#         else:
#             last.Lk['up'] += delta_up
#             last.sigma['up'] = _sigma(last, 'up')

#     # IV) atividade downstream (se estiver elongando)
#     first = RNAP_list[0]
#     if first.t_elongating:
#         L_down = first.Lk0['down'] * modelP.dna.n
#         D_dn = DLk_TopoI(L_down, first, modelP, loc="down")
#         G_dn = DLk_Gyrase(L_down, first, modelP, loc="down")
#         G_sp_dn = DLk_Gyrase("spec", first, modelP, loc="down")

#         delta_down = D_dn + G_dn + G_sp_dn
#         if delta_down != 0:
#             modelP.gene.Lk_domain += delta_down
#             first.Lk['down'] += delta_down
#             first.sigma['down'] = _sigma(first, 'down')


# def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory): ESSE AQUI ESTÁ QUASE, MAS DEMORA UM POUCO PARA INICIALIZAR A TOPO1
#     """promoter escape => RNAP is now in elongating mode"""

#     RNAP_list[-1].t_elongating = True
#     RNAP_list[-1].tesc = traj.time
#     RNAP_list[-1].topo1_bound = True  # NOVA: TopoI inicialmente ligada após escape

#     traj.esc_times["mean"] = (
#         traj.esc_times["n"] * traj.esc_times["mean"]
#         + RNAP_list[-1].tesc
#         - RNAP_list[-1].tocf
#     ) / (traj.esc_times["n"] + 1)
#     traj.esc_times["n"] += 1

#     traj.initiation_times["mean"] = (
#         traj.initiation_times["n"] * traj.initiation_times["mean"]
#         + RNAP_list[-1].tesc
#         - RNAP_list[-1].tb
#     ) / (traj.initiation_times["n"] + 1)
#     traj.initiation_times["n"] += 1

#     # UPDATING DOWNSTREAM RNAP
#     # change Lk because the new elongating RNA becomes a barrier
#     # one can actually check the conservation of the Lk
#     if len(RNAP_list) > 1:
#         RNAP_list[-2].Lk0['up'] -= RNAP_list[-1].Lk0['up']
#         RNAP_list[-2].Lk['up'] = (1 + RNAP_list[-2].sigma['up']) * RNAP_list[-2].Lk0['up']

#     return


# # 2. NOVA FUNÇÃO: topo1_unbinding_stage()
# def topo1_unbinding_stage(modelP: ModelParam, RNAP_list):
#     """TopoI unbinding from elongating RNAPs"""
    
#     # Se k_unbinding não existe ou é 0, mantém comportamento original
#     if not hasattr(modelP, 'k_unbind') or modelP.k_unbind == 0:
#         return
    
#     # Processa cada RNAP para possível unbinding
#     for rnap in RNAP_list:
#         if rnap.t_elongating and hasattr(rnap, 'topo1_bound') and rnap.topo1_bound:
#             # Probabilidade de unbinding por iteração
#             p_unbinding = modelP.k_unbind * modelP.coarse_g.tau_0
#             if np.random.uniform() < p_unbinding:
#                 rnap.topo1_bound = False
    
#     return


# # 3. MODIFICAÇÃO NA FUNÇÃO topo_stage_RNAPpresent()
# def topo_stage_RNAPpresent(modelP: ModelParam, RNAP_list):
#     """TopoI and gyrase activity in the presence of at least one DNA-bound RNAP"""

#     # UPSTREAM
#     # non-specific activities
#     DTopoI, DGyrase = 0, 0
#     if not RNAP_list[-1].t_elongating:
#         # the most upstream RNAP (at the promoter) is not a barrier
#         if len(RNAP_list) == 1:
#             # topoisomerases can act anywhere along the domain
#             domain_length_topo = modelP.gene.L_domain
#             DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#             DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
#         else:
#             # the second RNAP is a barrier and we consider activity upstream
#             domain_length_topo = RNAP_list[-2].Lk0['up'] * modelP.dna.n
#             DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#             DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
#     else:
#         # the most upstream RNAP is a barrier and we consider activity upstream
#         domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
#         DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
#         DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)

#     # MODIFICAÇÃO: specific activity in the presence of transcription (only TopoI)
#     DTopoI_spec = 0
#     if len(RNAP_list) > 1 or RNAP_list[-1].t_elongating:
#         # Verifica se TopoI está ligada à RNAP antes de aplicar atividade específica
#         if (hasattr(RNAP_list[-1], 'topo1_bound') and RNAP_list[-1].topo1_bound) or \
#            not hasattr(RNAP_list[-1], 'topo1_bound'):  # Compatibilidade com código antigo
#             DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)

#     if DTopoI != 0 or DGyrase != 0 or DTopoI_spec != 0:  # updating topo properties
#         modelP.gene.Lk_domain += DTopoI + DGyrase + DTopoI_spec

#         if not RNAP_list[-1].t_elongating:
#             # properties of non-elongating RNAP are dictated by its downstream RNAP (if it exists)
#             if len(RNAP_list) == 1:
#                 RNAP_list[0].sigma['up'] = (
#                     modelP.gene.Lk_domain - modelP.gene.Lk0_domain
#                 ) / modelP.gene.Lk0_domain
#             else:
#                 RNAP_list[-2].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
#                 RNAP_list[-2].sigma['up'] = _sigma(RNAP_list[-2], 'up')
#                 RNAP_list[-1].sigma['up'] = RNAP_list[-2].sigma['up']

#             RNAP_list[-1].Lk['up'] = (1 + RNAP_list[-1].sigma['up']) * RNAP_list[-1].Lk0['up']
#             RNAP_list[-1].sigma['down'] = RNAP_list[-1].sigma['up']
#             RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[-1].Lk0['down']
#         else:
#             # RNAP is a barrier
#             RNAP_list[-1].Lk['up'] += DTopoI + DGyrase + DTopoI_spec
#             RNAP_list[-1].sigma['up'] = _sigma(RNAP_list[-1], 'up')

#     # DOWNSTREAM
#     DTopoI_down, DGyrase_down, DGyrase_spec = 0, 0, 0
#     if RNAP_list[0].t_elongating:
#         # if non elongating, this means a single non-elongating RNAP => treated at the upstream level
#         domain_length_topo = RNAP_list[0].Lk0['down'] * modelP.dna.n
#         DTopoI_down = DLk_TopoI(domain_length_topo, RNAP_list[0], modelP, loc="down")
#         DGyrase_down = DLk_Gyrase(domain_length_topo, RNAP_list[0], modelP, loc="down")
#         DGyrase_spec = DLk_Gyrase("spec", RNAP_list[0], modelP, loc="down")

#         modelP.gene.Lk_domain += DTopoI_down + DGyrase_down + DGyrase_spec
#         RNAP_list[0].Lk['down'] += DTopoI_down + DGyrase_down + DGyrase_spec
#         RNAP_list[0].sigma['down'] = _sigma(RNAP_list[0], 'down')

#     return

def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
    """promoter escape => RNAP is now in elongating mode"""

    RNAP_list[-1].t_elongating = True
    RNAP_list[-1].tesc = traj.time

    #rnap = RNAP_list[-1]
    #rnap.topo_spec.bind()


    # MODIFICAÇÃO: Topo1 sempre presente desde o início (comportamento do modelo antigo)
    #RNAP_list[-1].topo_spec.bound = True 
    # se for comportamento antigo, força topo1, senão não mexe aqui:

    # if not modelP.enable_topo1_unbinding: tirei aqui 18/07/16:00
    #     RNAP_list[-1].topo1_bound = True

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


# def topo1_unbinding_stage(modelP: ModelParam, RNAP_list):
#     """TopoI unbinding from elongating RNAPs"""
    
#     # MODIFICAÇÃO: Desabilita unbinding por padrão para manter comportamento original
#     # Se k_unbinding não existe, é 0, ou se queremos comportamento original
#     if (not hasattr(modelP, 'k_unbind') or 
#         modelP.k_unbind == 0 or 
#         not hasattr(modelP, 'enable_topo1_unbinding') or 
#         not modelP.enable_topo1_unbinding):
#         return
    
#     # Processa cada RNAP para possível unbinding (só se explicitamente habilitado)
#     for rnap in RNAP_list:
#         if rnap.t_elongating and hasattr(rnap, 'topo1_bound') and rnap.topo1_bound:
#             # Probabilidade de unbinding por iteração
#             p_unbinding = modelP.k_unbind * modelP.coarse_g.tau_0
#             if np.random.uniform() < p_unbinding:
#                 rnap.topo1_bound = False
    
#     return

def topo1_unbinding_stage(modelP: ModelParam, RNAP_list):
    """TopoI unbinding from elongating RNAPs"""
    #print(f"[DEBUG] topo1_unbinding_stage called: k_unbind={modelP.k_unbind}, "
    #f"enable={modelP.enable_topo1_unbinding}, tau_0={modelP.coarse_g.tau_0}")

    # if (not hasattr(modelP, 'k_unbind') or
    #     modelP.k_unbind == 0 or
    #     not hasattr(modelP, 'enable_topo1_unbinding') or
    #     not modelP.enable_topo1_unbinding):
    #     #print("[DEBUG] unbinding guard tripped, exiting")
    #     return
    for rnap in RNAP_list:
        #if rnap.t_elongating and rnap.topo1_bound:
        if rnap.topo_spec.bound:
            p_unbinding = modelP.k_unbind * modelP.coarse_g.tau_0
            #print(f"[DEBUG] testing RNAP[{id(rnap)}]: p_unbind={p_unbinding:.3e}")
            if p_unbinding > 0 and np.random.uniform() < p_unbinding:
                #print(f"[DEBUG]   -> unbound topo1 from RNAP[{id(rnap)}]")
                rnap.topo_spec.unbind()
                _stats['unbind_events'] += 1



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
            DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)
    else:
        # the most upstream RNAP is a barrier and we consider activity upstream
        domain_length_topo = RNAP_list[-1].Lk0['up'] * modelP.dna.n
        DTopoI = DLk_TopoI(domain_length_topo, RNAP_list[-1], modelP)
        DGyrase = DLk_Gyrase(domain_length_topo, RNAP_list[-1], modelP)

    # MODIFICAÇÃO: specific activity - assume topo1 sempre presente por padrão
    DTopoI_spec = 0
    if len(RNAP_list) > 1 or RNAP_list[-1].t_elongating:
        # Verifica se TopoI está ligada à RNAP
        # Por padrão, assume que está sempre ligada (comportamento modelo antigo)
        # topo1_is_bound = True
        # #if hasattr(RNAP_list[-1], 'topo1_bound'):
        # topo1_is_bound = RNAP_list[-1].topo1_bound
        
        # if topo1_is_bound:
        #     DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)
        if RNAP_list[-1].topo_spec.bound:
            DTopoI_spec = DLk_TopoI("spec", RNAP_list[-1], modelP)

    _stats['spec_events'] += (DTopoI_spec if isinstance(DTopoI_spec, int) else int(DTopoI_spec))

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
            RNAP_list[-1].sigma['down'] = RNAP_list[-1].sigma['up']
            RNAP_list[-1].Lk['down'] = (1 + RNAP_list[-1].sigma['down']) * RNAP_list[-1].Lk0['down']
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


# FUNÇÃO ADICIONAL: Para inicializar RNAP com topo1 sempre presente
#def initialize_RNAP_with_topo1(RNAP_instance):
    """Inicializa uma instância de RNAP com topo1 sempre ligada (comportamento modelo antigo)"""
    RNAP_instance.topo1_bound = True
    return RNAP_instance


# FUNÇÃO ADICIONAL: Para controlar comportamento do modelo
#def set_model_behavior(modelP: ModelParam, use_original_behavior=True):
    """
    Define o comportamento do modelo:
    - use_original_behavior=True: Topo1 sempre presente (modelo antigo)
    - use_original_behavior=False: Permite unbinding de topo1 (modelo novo)
    """
    if use_original_behavior:
        modelP.enable_topo1_unbinding = False
        modelP.k_unbind = 0
    else:
        modelP.enable_topo1_unbinding = True
        # k_unbind deve ser definido conforme necessário
    
    return modelP

