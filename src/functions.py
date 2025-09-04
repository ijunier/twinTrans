import numpy as np, os, sys
import math
from param_var import *

"""
    Functions used by ../bin/twin.py    
"""
# _stats = {
#     'bind_oc': 0,
#     'bind_bind': 0,
#     'spec_events': 0,
#     'unbind_events': 0,
#     'iter': 0
# }

# def _dump_stats(label):
#     print(f"\n=== STATS {label} ===")
#     for k,v in _stats.items():
#         print(f"{k}: {v}")

def write_elongation_traj(traj: Trajectory, RNAP_list, simuP: SimuParam, header=False):
    """
    Write time series of number of elongating RNAPs.
    File: <output_folder>/traj_elongating.txt
    Columns: time \t N_elong
    """
    fi = simuP.fo_out + "/traj_elongating.txt"
    if header:
        with open(fi, "w") as out:
            out.write("time\tN_elong\n")
        return
    with open(fi, "a") as out:
        N_elong = sum(1 for r in RNAP_list if r.t_elongating)
        out.write(f"{traj.time}\t{N_elong}\n")
    return

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
        
        # New: TopoI unbinding
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

    # --- sampling for elongation counts ---
    traj.elong_counts = []  # store (time, N_elong)
    sample_every_secs = 1.0
    sample_every_iters = max(1, int(sample_every_secs / modelP.coarse_g.tau_0))

    write_transcripts_on_the_fly(traj, simuP, header=True)
    while traj.niter < simuP.Niterations and traj.Ntranscripts < simuP.Ntranscripts_max:
        verbosing(simuP, traj, RNAP_list)
        #update_presence_topo1_spec(modelP, RNAP_list)

        # Binding
        if traj.niter == nextevent2iter["b"]:
            
            binding_stage(modelP, RNAP_list, traj, nextevent2iter)

        # New: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # OC formation
        if traj.niter == nextevent2iter["oc"]:
            oc_formation_stage(modelP, RNAP_list, traj, nextevent2iter)
            

        # New: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        
        # Promoter escape
        if traj.niter == nextevent2iter["esc"]:
            escape_stage(modelP, RNAP_list, traj)

        # New: TopoI unbinding
        topo1_unbinding_stage(modelP, RNAP_list)


        # Topoisomerases
        topo_stage(modelP, RNAP_list)

        # New: TopoI unbinding
        #topo1_unbinding_stage(modelP, RNAP_list)

        # Elongation
        elongation_stage(modelP, RNAP_list)

        # Termination
        termination_stage(modelP, simuP, RNAP_list, traj)

        # if traj.time >= modelP.promoter._t_off and not RNAP_list:
        #     break

        # --- sample elongating RNAPs at configured frequency ---
        if not traj.niter % sample_every_iters:
            N_elong = sum(1 for r in RNAP_list if r.t_elongating)
            traj.elong_counts.append((traj.time, N_elong))
            write_elongation_traj(traj, RNAP_list, simuP)


        traj.niter += 1
        traj.time = traj.niter * modelP.coarse_g.tau_0


    #_dump_stats("multiple_transcripts")
    return #traj


# Transcription stages
def binding_stage(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter):
    """Binding of the RNAP at the promoter"""

    # if traj.time < modelP.promoter._t_off:
    #     kb_eff = modelP.promoter.kb
    # else:
    #     kb_eff = 0.0  novo




    # if kb_eff > 0:
    # # escala = 1/λ para experimento exponencial
    #     Z = np.random.exponential(scale=1.0/kb_eff)
    #     next_b = traj.niter + int(Z / modelP.coarse_g.tau_0)
    # else:
    #     # após shutdown, nunca mais agenda binding
    #     next_b = sys.maxsize
        
    # nextevent2iter["b"] = max(traj.niter + 1, next_b)  novo

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
        #rnap.topo_spec.bind()

        #rnap.is_on = traj.time < modelP.promoter._t_off novo
        
        # _stats['bind_bind'] += 1 #DEBUG Counting how many binding events of topo1 occured
        
        
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
    # Z = np.random.exponential(scale=1.0/kb_eff, size=None) if kb_eff>0 else float('inf')
    # nextevent2iter["b"] = np.max(
    #     (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
    # )  # at least niter + 1 NOVO

          # NEXT BINDING TRIAL (INDEPENDENT WHETHER BINDING HAS OCCURED OR NOT)
    Z = np.random.exponential(scale=modelP.promoter.kb_s, size=None)
    nextevent2iter["b"] = np.max(
        (traj.niter + 1, traj.niter + int(Z / modelP.coarse_g.tau_0))
    )  # at least niter + 1


    return




def oc_formation_stage(modelP: ModelParam, RNAP_list, traj: Trajectory, nextevent2iter): #MODELO NOVO
    """OC formation if sigma <= threshold"""
       

    if RNAP_list[-1].sigma['up'] <= modelP.promoter.sigma_o:
        #RNAP_list[-1].topo_spec.bind()
        Topo1Spec.bind()
        
        # _stats['bind_oc'] += 1  #DEBUG: Counting how many events of binding in OC Formation happened

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



def topo_stage(modelP: ModelParam, RNAP_list):
    """TopoI and gyrase activity"""

    if RNAP_list:
        topo_stage_RNAPpresent(modelP, RNAP_list)
    else:
        topo_stage_RNAPabsent(modelP)

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


# def termination_stage(modelP: ModelParam, simuP: SimuParam, RNAP_list, traj: Trajectory):
#     """Termination stage: transcript production by the most downstream RNAP"""

#     if RNAP_list and RNAP_list[0].X >= modelP.gene.term:
#         #traj.Ntranscripts += 1

#         # 1) compute speed of this RNAP
#         r = RNAP_list[0]




#         # --- 1) calcule tempo de elongação desde o escape
#         # elong_time = traj.time - r.tesc
#         # # --- 2) velocidade em nt/s (gene.L é o comprimento em bp)
#         # speed = modelP.gene.L / elong_time

#         # # --- 3) classifique on/off pelo instante de término
#         # if traj.time <= modelP.promoter._t_off:
#         #     traj.vel_on .append(speed)
#         #     traj.n_before_off += 1
#         # else:
#         #     traj.vel_off.append(speed)
#         #     traj.n_after_off  += 1

#         # if traj.time < modelP.promoter._t_off:
#         #     traj.n_before_off += 1
#         # else:
#         #     traj.n_after_off  += 1

#         traj.termination_times.append(traj.time)
#         traj.termination_escape_times.append(r.tesc)




#         # 3) increment transcript count
#         traj.Ntranscripts += 1


#         # ELONGATION TIME
#         traj.elongation_times["mean"] = (
#             traj.elongation_times["n"] * traj.elongation_times["mean"]
#             + traj.time
#             - RNAP_list[0].tesc
#         ) / (traj.elongation_times["n"] + 1)
#         traj.elongation_times["n"] += 1

#         # PRODUCTION TIME
#         if traj.time_last_prod > 0:
#             traj.prod_times["mean"] = (
#                 traj.prod_times["n"] * traj.prod_times["mean"]
#                 + traj.time
#                 - traj.time_last_prod
#             ) / (traj.prod_times["n"] + 1)
#             traj.prod_times["n"] += 1
#         traj.time_last_prod = traj.time

#         # TRANSCRIPT TERMINATION
#         # 1. We update the upstream RNAP, if it exists        
#         if len(RNAP_list) > 1: # at least 2 RNAPs: upating of the upstream RNAP (index = 1)
#             RNAP_list[1].Lk0['down'] += RNAP_list[0].Lk0['down']
#             if RNAP_list[1].t_elongating:
#                 # only downstream properties are updated                
#                 RNAP_list[1].Lk['down'] += RNAP_list[0].Lk['down']
#                 RNAP_list[1].sigma['down'] = _sigma(RNAP_list[1], 'down')
#             else:
#                 # both upstrean and dowsntream properties are updtaed from sigma of the domain
#                 RNAP_list[1].sigma['up'] = (
#                     modelP.gene.Lk_domain - modelP.gene.Lk0_domain
#                 ) / modelP.gene.Lk0_domain
#                 RNAP_list[1].Lk['up'] = (1 + RNAP_list[1].sigma['up']) * RNAP_list[1].Lk0['up']

#                 RNAP_list[1].sigma['down'] = RNAP_list[1].sigma['up']
#                 RNAP_list[1].Lk['down'] = (1 + RNAP_list[1].sigma['down']) * RNAP_list[1].Lk0['down']
#         # 2: We remove the RNAP
#         del RNAP_list[0]

#         #if not traj.Ntranscripts % simuP.Nevery_transcripts: (original)
#         if traj.Ntranscripts == 1 or not traj.Ntranscripts % simuP.Nevery_transcripts:
#             write_transcripts_on_the_fly(traj, simuP)

        

#     return

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



def topo1_unbinding_stage(modelP, RNAP_list):
    """Unbinding of topo1
    """
    # if there is TOPO1, it is going to unbind topo1
    if Topo1Spec.is_bound:
        # probability of unbiding each interaction
        p_unbind = modelP.k_unbind * modelP.coarse_g.tau_0
        if p_unbind > 0 and np.random.uniform() < p_unbind:
            Topo1Spec.unbind()           # Unbind TopoI 
            #_stats['unbind_events'] += 1  #counting how many unbind events



def topo_stage_RNAPpresent(modelP: ModelParam, RNAP_list):
    """TopoI and gyrase activity in the presence of at least one DNA-bound RNAP"""
    
    # UPSTREAM
    # non-specific activitiesdef+
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

    # MODIFICATION
    DTopoI_spec = 0

    rnap_up = RNAP_list[-1]
    if Topo1Spec.is_bound():  
        DTopoI_spec = DLk_TopoI("spec", rnap_up, modelP)
    else:
        DTopoI_spec = 0


    #_stats['spec_events'] += (DTopoI_spec if isinstance(DTopoI_spec, int) else int(DTopoI_spec))  # DEBUG:Couting how many spec events

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




def escape_stage(modelP: ModelParam, RNAP_list, traj: Trajectory):
    """promoter escape => RNAP is now in elongating mode"""
    #rnap = RNAP_list[-1]
    RNAP_list[-1].t_elongating = True
    RNAP_list[-1].tesc = traj.time
    #rnap.topo_spec.bind()

    #RNAP_list[-1].topo_spec.bound = True #new 18/07

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












