import numpy as np
from scipy import spatial as sp
#import matplotlib.pyplot as plt
import random
import subprocess
import sys
import time





def quat_mul(qa, qb):
    r0 = qa[0] * qb[0] - qa[1] * qb[1] - qa[2] * qb[2] - qa[3] * qb[3]
    r1 = qa[0] * qb[1] + qa[1] * qb[0] + qa[2] * qb[3] - qa[3] * qb[2]
    r2 = qa[0] * qb[2] - qa[1] * qb[3] + qa[2] * qb[0] + qa[3] * qb[1]
    r3 = qa[0] * qb[3] + qa[1] * qb[2] - qa[2] * qb[1] + qa[3] * qb[0]
    return np.array([r0, r1, r2, r3])

def quat_axis_angle(axis, angle):
    qw = np.cos(0.5 * angle)
    qx = axis[0] * np.sin(0.5 * angle)
    qy = axis[1] * np.sin(0.5 * angle)
    qz = axis[2] * np.sin(0.5 * angle)
    return np.array([qw, qx, qy, qz])

def mag(arg):
    """
    magnitude of 3d vector
    """
    return np.sqrt(arg[0] * arg[0] + arg[1] * arg[1] + arg[2] * arg[2])



def unit_vec(vector):
    """
    returns vector converted to a unit vector
    """
    m = mag(vector)
    return np.array([vector[0] / m, vector[1] / m, vector[2] / m])

def quat_norm(q):
    m = np.linalg.norm(q)
    return np.array([q[0] / m, q[1] / m, q[2] / m, q[3] / m])

def rotation(xin, axis, angle):
    """
    y = rotate "x" by "angle" about "axis"  axis needs to be unit vector
    uses Rodrigue's Rotation Formula:
    y = x cos(t) + (k cross x)sin(t) + k(k.x)(1-cos(t))
    where t = angle, k = axis

    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    """
    axis = unit_vec(axis)
    kcrossx = np.cross(axis, xin)
    kdotx = np.dot(axis, xin)

    # print kcrossx
    # print kdotx
    yout = np.array([0.0, 0.0, 0.0])

    yout[0] = xin[0] * np.cos(angle) + kcrossx[0] * np.sin(angle) \
              + axis[0] * kdotx * (1.0 - np.cos(angle))
    yout[1] = xin[1] * np.cos(angle) + kcrossx[1] * np.sin(angle) \
              + axis[1] * kdotx * (1.0 - np.cos(angle))
    yout[2] = xin[2] * np.cos(angle) + kcrossx[2] * np.sin(angle) \
              + axis[2] * kdotx * (1.0 - np.cos(angle))
    return yout


# def get_min_distance_mean_all(name):
#     """ Get the min distance between P and CA beads in the pdb file"""
#     pdb_file = list(open(name, "r"))
#
#     aa_dict = {
#         "ALA": list(),
#         "ARG": list(),
#         "ASN": list(),
#         "ASP": list(),
#         "CYS": list(),
#         "GLN": list(),
#         "GLU": list(),
#         "GLY": list(),
#         "HIS": list(),
#         "ILE": list(),
#         "LEU": list(),
#         "LYS": list(),
#         "MET": list(),
#         "PHE": list(),
#         "PRO": list(),
#         "SER": list(),
#         "THR": list(),
#         "TRP": list(),
#         "TYR": list(),
#         "VAL": list(),
#     }
#
#     Ps = list()
#     old_rnum = 1
#     for line in pdb_file:
#         split_line = line.split()
#
#         if split_line[0] == "ATOM":
#             # get amino acid, we know chains J and K are DNA
#             if split_line[-2] != "J" and split_line[-2] != "K":
#                 rnum = split_line[5]
#                 aa = split_line[3]
#
#                 if old_rnum != rnum:
#                     aa_dict[aa].append(list())
#
#                 chain = split_line[-2]
#
#                 x = split_line[-7]
#                 y = split_line[-6]
#                 z = split_line[-5]
#
#                 bead = np.array([x, y, z], dtype=float)
#
#                 aa_dict[aa][-1].append(bead)
#                 old_rnum = rnum
#
#             # get dna phosphates
#             if split_line[2] == "P":
#                 aa = split_line[3]
#
#                 chain = split_line[-2]
#
#                 x = split_line[-7]
#                 y = split_line[-6]
#                 z = split_line[-5]
#
#                 bead = np.array([x, y, z], dtype=float)
#
#                 Ps.append(bead)
#
#     # search through all P-Ca pairs
#     # put into pure numpy arrays to make it quicker
#     Ps = np.array(Ps)
#
#     # do the same for the ref frames
#     dna_data = list(open("ref_frames.dat"))
#     n_dna = int(dna_data[0].split()[0])
#
#     dna_beads = list()
#     line = 1
#     for n in range(n_dna):
#         line = line + 1
#         coods = dna_data[line].split()
#         r = np.array(
#             [str(coods[0]),
#              str(coods[1]),
#              str(coods[2])],
#             dtype=float)
#         line = line + 1
#
#         line = line + 1
#
#         line = line + 1
#
#         line = line + 1
#
#         dna_beads.append(r)
#
#     dna_beads = np.array(dna_beads)
#
#     # loop over all amino acids
#     min_Ps = list()
#     min_DNAs = list()
#     for key in aa_dict:
#         value = aa_dict[key]
#         # check it is not empty
#         if len(value) > 1:
#             # average each component
#             values = list()
#
#             for val in value:
#                 valmeans = np.mean(np.array(val), axis=0)
#                 values.append(valmeans)
#
#             dist_arr = sp.distance.cdist(Ps, np.array(values))
#
#             min_d = np.min(dist_arr)
#             min_Ps.append(min_d)
#             print("P  and ", key, " = ", min_d)
#
#             dist_arr = sp.distance.cdist(dna_beads, np.array(values))
#             min_d = np.min(dist_arr)
#             min_DNAs.append(min_d)
#             print("bp and ", key, " = ", min_d)
#
#     mean_p = np.mean(min_Ps)
#     mean_DNA = np.mean(min_DNAs)
#
#     print(mean_p)
#     print(mean_DNA)
#
#
# def get_min_distance_CAs(name):
#     """ Get the min distance between P and CA beads in the pdb file"""
#     pdb_file = list(open(name, "r"))
#
#     aa_dict = {
#         "ALA": list(),
#         "ARG": list(),
#         "ASN": list(),
#         "ASP": list(),
#         "CYS": list(),
#         "GLN": list(),
#         "GLU": list(),
#         "GLY": list(),
#         "HIS": list(),
#         "ILE": list(),
#         "LEU": list(),
#         "LYS": list(),
#         "MET": list(),
#         "PHE": list(),
#         "PRO": list(),
#         "SER": list(),
#         "THR": list(),
#         "TRP": list(),
#         "TYR": list(),
#         "VAL": list(),
#     }
#
#     Ps = list()
#     for line in pdb_file:
#         split_line = line.split()
#
#         if split_line[0] == "ATOM":
#             # get amino acid, we know chains J and K are DNA
#             if split_line[2] == "CA":
#                 aa = split_line[3]
#
#                 chain = split_line[-2]
#
#                 x = split_line[-7]
#                 y = split_line[-6]
#                 z = split_line[-5]
#
#                 bead = np.array([x, y, z], dtype=float)
#
#                 aa_dict[aa].append(bead)
#
#             # get dna phosphates
#             if split_line[2] == "P":
#                 aa = split_line[3]
#
#                 chain = split_line[-2]
#
#                 x = split_line[-7]
#                 y = split_line[-6]
#                 z = split_line[-5]
#
#                 bead = np.array([x, y, z], dtype=float)
#
#                 Ps.append(bead)
#
#     # search through all P-Ca pairs
#     # put into pure numpy arrays to make it quicker
#     Ps = np.array(Ps)
#
#     # do the same for the ref frames
#     dna_data = list(open("ref_frames.dat"))
#     n_dna = int(dna_data[0].split()[0])
#
#     dna_beads = list()
#     line = 1
#     for n in range(n_dna):
#         line = line + 1
#         coods = dna_data[line].split()
#         r = np.array(
#             [str(coods[0]),
#              str(coods[1]),
#              str(coods[2])],
#             dtype=float)
#         line = line + 1
#
#         line = line + 1
#
#         line = line + 1
#
#         line = line + 1
#
#         dna_beads.append(r)
#
#     dna_beads = np.array(dna_beads)
#
#     # loop over all amino acids
#     min_Ps = list()
#     min_DNAs = list()
#     for key in aa_dict:
#         values = aa_dict[key]
#         # check it is not empty
#         if len(values) > 1:
#             # average each component
#
#             dist_arr = sp.distance.cdist(Ps, np.array(values))
#
#             min_d = np.min(dist_arr)
#             min_Ps.append(min_d)
#             print("P  and ", key, " = ", min_d)
#
#             dist_arr = sp.distance.cdist(dna_beads, np.array(values))
#             min_d = np.min(dist_arr)
#             min_DNAs.append(min_d)
#             print("bp and ", key, " = ", min_d)
#
#     mean_p = np.mean(min_Ps)
#     mean_DNA = np.mean(min_DNAs)
#
#     print(mean_p)
#     print(mean_DNA)


def exyz_to_q(ex, ey, ez):
    """ taken from LAMMPS source code, converts ex,ey,ez frame vectors to quaternion orientations """

    # squares of quaternion components

    q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0)
    q1sq = q0sq - 0.5 * (ey[1] + ez[2])
    q2sq = q0sq - 0.5 * (ex[0] + ez[2])
    q3sq = q0sq - 0.5 * (ex[0] + ey[1])

    q = np.array([0.0, 0.0, 0.0, 0.0])
    # some component must be greater than 1/4 since they sum to 1
    # compute other components from it

    if q0sq >= 0.25:
        q[0] = np.sqrt(q0sq)
        q[1] = (ey[2] - ez[1]) / (4.0 * q[0])
        q[2] = (ez[0] - ex[2]) / (4.0 * q[0])
        q[3] = (ex[1] - ey[0]) / (4.0 * q[0])
    elif q1sq >= 0.25:
        q[1] = np.sqrt(q1sq)
        q[0] = (ey[2] - ez[1]) / (4.0 * q[1])
        q[2] = (ey[0] + ex[1]) / (4.0 * q[1])
        q[3] = (ex[2] + ez[0]) / (4.0 * q[1])
    elif q2sq >= 0.25:
        q[2] = np.sqrt(q2sq)
        q[0] = (ez[0] - ex[2]) / (4.0 * q[2])
        q[1] = (ey[0] + ex[1]) / (4.0 * q[2])
        q[3] = (ez[1] + ey[2]) / (4.0 * q[2])
    elif q3sq >= 0.25:
        q[3] = np.sqrt(q3sq)
        q[0] = (ex[1] - ey[0]) / (4.0 * q[3])
        q[1] = (ez[0] + ex[2]) / (4.0 * q[3])
        q[2] = (ez[1] + ey[2]) / (4.0 * q[3])

    norm = np.linalg.norm(q)
    q = q / norm

    return q


def q_to_ey(q):
    """ taken from LAMMPS source code converts q to ey space frame vector """

    ey = (2.0 * (q[1] * q[2] - q[0] * q[3]),
          q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
          2.0 * (q[2] * q[3] + q[0] * q[1]))

    return np.array(ey)


def q_to_exyz(q):
    ex = (q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3],
          2.0 * (q[1] * q[2] + q[0] * q[3]),
          2.0 * (q[1] * q[3] - q[0] * q[2]))

    ey = (2.0 * (q[1] * q[2] - q[0] * q[3]),
          q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3],
          2.0 * (q[2] * q[3] + q[0] * q[1]))

    ez = (2.0 * (q[1] * q[3] + q[0] * q[2]),
          2.0 * (q[2] * q[3] - q[0] * q[1]),
          q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3])

    return np.array(ex), np.array(ey), np.array(ez)

def lj_pot(r,s,e):
    return 4*e*((s/r)**12 - (s/r)**6)

def lk_KA_test_potential(r,s,e,l):

    if r < 2.0**(1.0/6.0)*s:
        return lj_pot(r,s,e) + (1-l)*e
    else:
        return l*lj_pot(r,s,e)



def get_single_nucl(WIDTH, P_ANGLE):
    # load in the pdb file

    pdb_file = list(open("clust.pdb", "r"))

    # parse the data in the file
    # make list of x,y,z, amino acid, protein backbone

    data = list()

    for line in pdb_file:
        split_line = line.split()

        if split_line[0] == "ATOM":
            if split_line[2] == "CA":
                aa = split_line[3]

                chain = split_line[-2]

                x = split_line[-7]
                y = split_line[-6]
                z = split_line[-5]

                bead = [aa, chain, np.array([x, y, z], dtype=float)]

                data.append(bead)

    # get the dna data

    dna_data = list(open("ref_frames.dat"))

    # parse the data file
    # first line:
    # xxx base-pairs
    #
    # each bp organised like:
    # ...     1 G-C   # ...1>J:1169_:[DG5]G - ...1>K:1590_:[DC3]C\n'
    #   152.0681   242.6043   116.5436  # origin\n'
    #    -0.5248     0.6823    -0.5090  # x-axis\n'
    #     0.8288     0.5460    -0.1226  # y-axis\n'
    #     0.1943    -0.4862    -0.8520  # z-axis\n'

    n_dna = int(dna_data[0].split()[0])

    dna_beads = list()
    line = 1
    for n in range(n_dna):
        bps = dna_data[line].split()[2]
        line = line + 1
        coods = dna_data[line].split()
        r = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        x = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        y = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        z = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1

        # need to convert to a quaternion
        q = exyz_to_q(x, y, z)

        dna_beads.append([bps, r, q])

    # now add in phosphates
    dna_full = list()

    for n in range(0, n_dna):
        dna_full.append(["D", *dna_beads[n]])
        # need to get the positions of the phosphates
        ex, ey, ez = q_to_exyz(dna_full[-1][-1])

        pos1 = (dna_full[-1][2] + WIDTH * rotation(ey, ez, P_ANGLE))

        pos2 = (dna_full[-1][2] - WIDTH * rotation(ey, ez, -P_ANGLE))

        dna_full.append(["P1", "_", pos1])
        dna_full.append(["P2", "_", pos2])

    # split into fixed and globular domains
    # first split into the 9 different proteins

    proteins = list()
    p = list()
    proteins.append(p)
    old_c = "A"
    for line in data:
        if line[1] != old_c:
            # new chain
            p_new = list()
            proteins.append(p_new)
            proteins[-1].append(line)

        else:
            # current chain
            proteins[-1].append(line)
        old_c = line[1]

    # we have a seperate list for each protein, need to split them into tail and core
    # use defenitions from doi: 10.1529/biophysj.106.083006
    # Histone, chain, bps
    # H3  A, E N 1–40
    # H4  B, F N 1–25
    # H2A C, G N 1–20
    # H2A C, G C 114–128
    # H2B D, H N 1–25

    # H1 is the first one in our stucture so increment all above letters by 1

    indexes = [[[1, 24], [95, 194]],
               [[1, 40]],
               [[1, 25]],
               [[1, 20], [114, 128]],
               [[1, 25]],
               [[1, 40]],
               [[1, 25]],
               [[1, 20], [114, 128]],
               [[1, 25]],
               ]

    alt_proteins = list()
    p = 0
    for protein in proteins:
        indx = indexes[p]

        new_p = list()

        # print(len(protein))

        if len(indx) == 1:  # just N tail
            Ntail = list()
            fix = list()
            for i in range(indx[0][0], indx[0][1] + 1):
                Ntail.append(protein[i - 1])

            for i in range(indx[0][1] + 1, len(protein) + 1):
                fix.append(protein[i - 1])

            new_p.append(Ntail)
            new_p.append(fix)

        else:  # N and C tail
            Ntail = list()
            fix = list()
            Ctail = list()
            for i in range(indx[0][0], indx[0][1] + 1):
                Ntail.append(protein[i - 1])

            for i in range(indx[0][1] + 1, indx[1][0]):
                fix.append(protein[i - 1])

            for i in range(indx[1][0], len(protein) + 1):
                Ctail.append(protein[i - 1])

            new_p.append(Ntail)
            new_p.append(fix)
            new_p.append(Ctail)

        alt_proteins.append(new_p)
        p = p + 1

    return alt_proteins, dna_full, n_dna


def add_nucl(alt_proteins, dna_full,initial_dna_len):
    """  Adds in another nucleosome
    :param alt_proteins:
    :param dna_full:
    :return: new alt_proteins, updated dna_full
    """

    # We have the current structure
    # Use the DNA beads to construct the required transformation

    # get the last dna_beads

    last_dna = dna_full[-3]
    last_dna_q = last_dna[3]
    last_dna_x = last_dna[2]

    first_dna = dna_full[0]
    first_dna_q = first_dna[3]
    first_dna_x = first_dna[2]

    # make new dna bead
    last_dna_ex, last_dna_ey, last_dna_ez = q_to_exyz(last_dna_q)

    new_dna_x = last_dna_x + np.array(last_dna_ez) * 3.4
    new_dna_q = quat_mul(quat_axis_angle(np.array(last_dna_ez), 30 * np.pi / 180.0), last_dna_q)

    # Create a transformation matrix that maps first dna to new dna
    # 4x4 transformation matrix
    # M2 = T M1
    M1 = np.zeros([4, 4])
    M2 = np.zeros([4, 4])

    M1[:3, :3] = np.array(q_to_exyz(first_dna_q)).transpose()
    M1[:3, 3] = first_dna_x
    M1[3, 3] = 1.0

    M2[:3, :3] = np.array(q_to_exyz(new_dna_q)).transpose()
    M2[:3, 3] = new_dna_x
    M2[3, 3] = 1.0

    # T = M2 * M1^-1
    T = np.matmul(M2, np.linalg.inv(M1))

    # Use T to transform all the beads beads in the current structure to the new structure
    new_alt_proteins = list()
    new_dna_full = list()

    # for protein in alt_proteins:
    #     for subsection in protein:
    #         for atom in subsection:
    #             M = np.zeros([4, 4])

    new_dna = list()


    for n in range(initial_dna_len):
        dna = dna_full[n]

        if dna[0] == "D":
            M = np.zeros([4, 4])
            M[:3, :3] = np.array(q_to_exyz(dna[3])).transpose()
            M[:3, 3] = dna[2]
            M[3, 3] = 1.0
            new_M = np.matmul(T, M)
            new_q = exyz_to_q(new_M[:3, 0], new_M[:3, 1], new_M[:3, 2])
            new_x = new_M[:3, 3]
            new_dna.append([dna[0], dna[1], new_x, new_q])
        else:
            M = np.zeros([4, 4])
            M[:3, 3] = dna[2]
            M[3, 3] = 1.0
            new_M = np.matmul(T, M)
            new_x = new_M[:3, 3]
            new_dna.append([dna[0], dna[1], new_x])
    dna_full = dna_full + new_dna

    # Make new proteins
    new_alt_proteins = list()
    for protein in alt_proteins:
        new_protein = list()
        for subsection in protein:
            new_subsection = list()
            for atom in subsection:
                M = np.zeros([4, 4])
                x = atom[2]
                M[:3, 3] = x
                M[3, 3] = 1.0
                new_M = np.matmul(T, M)
                new_x = new_M[:3, 3]
                new_atom = [atom[0], atom[1], new_x]
                new_subsection.append(new_atom)
            new_protein.append(new_subsection)
        new_alt_proteins.append(new_protein)

    return new_alt_proteins, dna_full

def add_nrl(dna_full, left_l,right_l,WIDTH,P_ANGLE):
    RISE=3.4
    TWIST=35*180/np.pi
    current_nrl=len(dna_full)/3
    print("new nrl" , current_nrl+left_l+right_l)
    print("current nrl", current_nrl)
    #diff = new_nrl-current_nrl

    #print("diff",diff)


    #end1 = int(np.floor(diff/2))
    #end2 = int(np.ceil(diff/2))

    #print(end1)
    #print(end2)


    # end to the start

    if left_l < 0:

        dna_full = dna_full[(-left_l)*3:]
        #print(dna_full)


    else:
        first_dna = dna_full[0]
        x = first_dna[2]
        q = first_dna[3]

        new_dna=list()

        for i in range(0,left_l):

            ex,ey,ez = q_to_exyz(q)

            new_x = x - ez*RISE
            new_q = quat_norm( quat_mul( quat_norm( quat_axis_angle(ez,-TWIST)),q))

            # need to get the positions of the phosphates
            ex, ey, ez = q_to_exyz(new_q)
            unit_vec(ey)
            unit_vec(ez)

            pos1 = (new_x - WIDTH * rotation(ey, ez, -P_ANGLE))

            pos2 = (new_x + WIDTH * rotation(ey, ez, P_ANGLE))

            new_dna.append(["P1", "_", pos1])
            new_dna.append(["P2", "_", pos2])

            new_dna.append(["D", "A-T",new_x, new_q])

            q = quat_norm( new_dna[-1][3])
            x = new_dna[-1][2]


        new_dna.reverse()

        dna_full = new_dna + dna_full



    if right_l < 0:

        dna_full = dna_full[:(right_l*3)]

        
        #print(dna_full)
    else:

        last_dna = dna_full[-3]
        x = last_dna[2]
        q = quat_norm(last_dna[3])

        for i in range(0,right_l):

            ex,ey,ez = q_to_exyz(q)

            unit_vec(ez)

            new_x = x + ez*RISE
            new_q = quat_norm( quat_mul( quat_norm(quat_axis_angle(ez,TWIST)),q))

            unit_vec(ey)
            dna_full.append(["D", "A-T",new_x, new_q])
            # need to get the positions of the phosphates
            ex, ey, ez = q_to_exyz(quat_norm( dna_full[-1][-1]))

            pos1 = (dna_full[-1][2] + WIDTH * rotation(ey, ez, P_ANGLE))

            pos2 = (dna_full[-1][2] - WIDTH * rotation(ey, ez, -P_ANGLE))

            dna_full.append(["P1", "_", pos1])
            dna_full.append(["P2", "_", pos2])

            q = quat_norm(dna_full[-3][3])
            x = dna_full[-3][2]



    #print(dna_full)
    return dna_full




def create_nucl(EPSILON_P,EPSILON_D, HPS,NUM_NUCLS, version, LH, seq_type,RESTART,ADD_NRL,left_l, right_l):



    # 601 sequence:
    seq601 = "ACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAG"
    seqTTAGGG = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA"
    seqAT = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    RISE = 3.4
    TWIST = 36 * np.pi / 180.0

    # Charges
    charge_dict = {
        "ALA": 0,
        "ARG": 1,
        "ASN": 0,
        "ASP": -1,
        "CYS": 0,
        "GLN": 0,
        "GLU": -1,
        "GLY": 0,
        "HIS": 0.5,
        "ILE": 0,
        "LEU": 0,
        "LYS": 1,
        "MET": 0,
        "PHE": 0,
        "PRO": 0,
        "SER": 0,
        "THR": 0,
        "TRP": 0,
        "TYR": 0,
        "VAL": 0,
    }

    # VdW radii
    vdw_dict = {
        "ALA": 5.04,
        "ARG": 6.56,
        "ASN": 5.68,
        "ASP": 5.58,
        "CYS": 5.48,
        "CYD": 5.48,
        "GLN": 6.02,
        "GLU": 5.92,
        "GLY": 4.50,
        "HIS": 6.08,
        "ILE": 6.18,
        "LEU": 6.18,
        "LYS": 6.36,
        "MET": 6.18,
        "PHE": 6.36,
        "PRO": 5.56,
        "SER": 5.18,
        "THR": 5.62,
        "TRP": 6.78,
        "TYR": 6.46,
        "VAL": 5.86,
    }

    hps_lambdas = {
        "ILE": 0.972973,
        "VAL": 0.891892,
        "LEU": 0.972973,
        "PHE": 1.0,
        "CYS": 0.594595,
        "MET": 0.837838,
        "ALA": 0.72973,
        "GLY": 0.648649,
        "THR": 0.675676,
        "SER": 0.594595,
        "TRP": 0.945946,
        "TYR": 0.864865,
        "PRO": 1.0,
        "HIS": 0.513514,
        "GLU": 0.459459,
        "GLN": 0.513514,
        "ASP": 0.378378,
        "ASN": 0.432432,
        "LYS": 0.513514,
        "ARG": 0.0
    }

    # masses
    mass_dict = {
        "ALA": 71.08,
        "ARG": 156.20,
        "ASN": 114.10,
        "ASP": 115.10,
        "CYS": 103.10,
        "GLN": 128.10,
        "GLU": 129.10,
        "GLY": 57.05,
        "HIS": 137.10,
        "ILE": 113.20,
        "LEU": 113.20,
        "LYS": 128.20,
        "MET": 131.20,
        "PHE": 147.20,
        "PRO": 97.12,
        "SER": 87.08,
        "THR": 101.10,
        "TRP": 186.20,
        "TYR": 163.20,
        "VAL": 99.07,
    }

    D_a_epsilon = EPSILON_D
    P_a_epsilon = EPSILON_P

    P_lambda = 1.0
    D_lambda = 1.0


    # use elastic network model

    if version == 1 or version == 3 or version == 4 or version ==5:

        ENM = True
    else :
        ENM = False

    # cut-off distance for creating the elastic network
    EN_CUT = 7.5


    if HPS:
        D_type = 21
        P_type = 22

    else:
        D_type = 41
        P_type = 42


    DNA_vdw = 8
    P_vdw  = 4


    D_a_sigma_dict = {}

    P_a_sigma_dict = {}



    WIDTH = 8.5
    RX = 11.0
    RY = 20.0
    RZ = 3.5

    P_ANGLE = 20.0*np.pi/180.0

    Q_CUT = "$Q"

    if not HPS:
        # get the KH params dicts
        file_A = open("KH_A.txt","r")

        KH_A_dict = {}

        for line in file_A:
            split = line.split()
            KH_A_dict[split[0] + " " + split[1]] = float(split[2])


        file_D = open("KH_D.txt","r")

        KH_D_dict = {}

        for line in file_D:
            split = line.split()
            KH_D_dict[split[0] + " " + split[1]] = float(split[2])


    if version==2:
        in_file = open("in.nucl_v2","w")
    else:
        if not RESTART:
            in_file = open("in.nucl", "w")
        else:
            in_file = open("in.nucl_restart", "w")

    # VARIABLES
    if version == 2:
        print("variable fname index nucl_v2.txt", file=in_file)
        print("variable simname index nucl_v2", file=in_file)
    else:
        print("variable fname index nucl.txt", file=in_file)
        print("variable simname index nucl", file=in_file)

    print("", file=in_file)
    print("variable a equal 8.0", file=in_file)
    print("variable l equal 1.0/${a}", file=in_file)
    print("variable Q equal 3.5*${a}", file=in_file)
    print("", file=in_file)
    print("# Initialization", file=in_file)
    print("units real", file=in_file)
    print("atom_style   hybrid ellipsoid angle charge", file=in_file)
    print("", file=in_file)
    print("boundary s s s", file=in_file)
    print("", file=in_file)
    if  not RESTART:
        print("log log.${simname}.txt", file=in_file)
        print("read_data ${fname}", file=in_file)
    else:
        print("log log.${simname}.txt append", file=in_file)
        print("read_restart nucl.restart.*", file=in_file)

    print("", file=in_file)

    aa_keys = list(mass_dict)
    aa_keys.sort()

    # Masses
    # "set type    X mass  M"
    index_dict = {}
    index_dict_rev = {}
    a=1
    for aa in aa_keys:
        print("set type " + str(a) + " mass " + str(mass_dict[aa]), file=in_file)
        index_dict[a] = aa
        index_dict_rev[aa] = a
        a = a + 1

    if not HPS:
        for aa in aa_keys:
            print("set type " + str(a) + " mass " + str(mass_dict[aa]), file=in_file)
            index_dict[a] = aa
            a = a + 1


    print("set type " + str(D_type) + " mass 650.0", file=in_file)
    print("set type " + str(P_type) + " mass 0.000001", file=in_file)

    if version == 5 and LH:
        if HPS:
            print("set type " + str(23) + " mass 100.0", file=in_file)
        else:
            print("set type " + str(43) + " mass 100.0", file=in_file)


    #print("pair_style ljlambda 0.125 " + str(Q_CUT) + " " +str(Q_CUT), file=in_file)
    # print("pair_style ljlambda "+str(1.0/DEBYE_LENGTH) +" " + str(Q_CUT) + " " +str(Q_CUT), file=in_file)
    print("pair_style ljlambda $l 20.34 $Q", file=in_file)
    print("dielectric  80.0", file=in_file)
    print("", file=in_file)
    print("bond_style hybrid harmonic harmonic/DNA zero", file=in_file)
    print("", file=in_file)
    print("bond_coeff 1 harmonic 10.0 3.8", file=in_file)
    print("bond_coeff 2 harmonic/DNA 0.0 0.0", file=in_file)
    print("bond_coeff 3 zero", file=in_file)

    # get_min_distance_mean_all("clust.pdb")
    # get_min_distance_CAs("clust.pdb")

    # get the single nucleosome structure
    # DNA and histones
    alt_proteins, dna_full, n_dna = get_single_nucl(WIDTH, P_ANGLE)

    if ADD_NRL:
        dna_full = add_nrl(dna_full,left_l,right_l,WIDTH,P_ANGLE)

    dna_len_inital = len(dna_full)

    current_nrl =  dna_len_inital/3

    # replace the dna sequence
    if seq_type == 2:
        # only replace inner 147
        start_index = (current_nrl - 147)/2
        end_index = current_nrl - start_index

        print(start_index)
        print(end_index)

        for iiii in range(147):
            base = seq601[iiii]
            if base == "T":
                base2 = "A"
            elif base == "A":
                base2 = "T"
            elif base == "C":
                base2 = "G"
            else:
                base2 = "C"

            dna_idx = int((iiii + start_index)*3)
            dna_full[dna_idx][1] = base + "-" + base2

    if seq_type == 3:
        # only replace inner 147
        start_index = (current_nrl - 147)/2
        end_index = current_nrl - start_index

        print(start_index)
        print(end_index)

        for iiii in range(147):
            base = seqTTAGGG[iiii]
            if base == "T":
                base2 = "A"
            elif base == "A":
                base2 = "T"
            elif base == "C":
                base2 = "G"
            else:
                base2 = "C"

            dna_idx = int((iiii + start_index)*3)
            dna_full[dna_idx][1] = base + "-" + base2
    
    if seq_type == 4:
        # replace all
        for iiii in range(int(len(dna_full)/3)):
            dna_idx = int(iiii * 3)
            dna_full[dna_idx][1] = "A-T"

    if seq_type == "random":
        for iiii in range(int(len(dna_full)/3)):
            dna_idx = int(iiii * 3)
            dna_full[dna_idx][1] = np.random.choice(["A-T","T-A","C-G","G-C"])

    if not LH:
        alt_proteins = alt_proteins[1:]

    nucls = [alt_proteins]


    for N in range(1,NUM_NUCLS):

        # now we add in an extra nucleosomes
        new_alt_proteins, dna_full = add_nucl(alt_proteins, dna_full, dna_len_inital)

        nucls = nucls + [new_alt_proteins]

    # count number of DNA basepairs
    n_dna = 0
    for dna in dna_full:
        if dna[0] == "D":
            n_dna = n_dna + 1

    # sort out the bonds
    bonds = list()

    # counters for bond and bead index
    b = 1
    i = 1
    if ENM:
        # more complex process required here
        # need to make a network of all globular beads
        # need to identify their nearest neighbours

        # first assign tails as normal, and group the globs into seperate structure


        C_blobs = list()
        L_blobs = list()
        for alt_proteins in nucls:
            core_glob = list()
            H1_glob = list()
            p = 1
            for protein in alt_proteins:

                s = 1

                for subsection in protein:
                    if s == 1:
                        for atom in subsection:
                            bonds.append([b, 1, i, i + 1])
                            b = b + 1
                            i = i + 1
                    elif s == 2:
                        for atom in subsection:
                            if LH:
                                if p == 1:
                                    H1_glob.append([i, atom[2]])
                                else:
                                    core_glob.append([i, atom[2]])
                            else:
                                core_glob.append([i, atom[2]])

                            i = i + 1


                    elif s == 3:
                        for atom in subsection:
                            bonds.append([b, 1, i - 1, i])
                            b = b + 1
                            i = i + 1

                    s = s + 1

                p = p + 1
            C_blobs.append(core_glob)
            L_blobs.append(H1_glob)

        if version == 5 and LH:

            # get the com and add in dummy bead
            C_coms = list()
            L_coms = list()

            for blob in L_blobs:
                pos = list()
                for atom in blob:
                    pos.append(atom[1])

                pos = np.array(pos)
                com = np.mean(pos, axis=0)
                print(com)
                L_coms.append(com)
            for blob in C_blobs:
                pos = list()
                for atom in blob:
                    pos.append(atom[1])

                pos = np.array(pos)
                com = np.mean(pos, axis=0)
                print(com)
                C_coms.append(com)
            dummy_bond_ids = list()
            dummy_bond_dist = list()
            dummy_atom = list()
            for cc in range(len(C_blobs)):
                L_blobs[cc].append([i, L_coms[cc]])
                dummy_atom.append(["X", "X", L_coms[cc]])
                i = i + 1
                C_blobs[cc].append([i, C_coms[cc]])
                dummy_atom.append(["X", "X", C_coms[cc]])
                dummy_bond_ids.append([i, i - 1])
                dummy_bond_dist.append(np.linalg.norm(C_coms[cc] - L_coms[cc]))
                i = i + 1

        # print(core_glob)


        if version == 3: # need to add in most of the nucleosomeal dna

            #need to get the atom id of the one before first dna.


            dna_nucl_sections = list()
            dna_nucl_atom_ids = list()
            nrl = int(dna_len_inital / 3)

            mid_point = int(np.floor(nrl / 2))

            start_point = int(mid_point - 142 / 2)
            end_point = int(mid_point + 142 / 2)

            dna_blobs = list()

            print(i)

            for N in range(NUM_NUCLS):
                dna_nucl_sections.append([start_point * 3 + N * nrl * 3, end_point * 3 + N * nrl * 3 + 2])

                dna_nucl_atom_ids.append([dna_nucl_sections[-1][0] + i, dna_nucl_sections[-1][1] + i])

                dna_glob = list()
                for i_idx,d_idx in zip(range(dna_nucl_atom_ids[-1][0], dna_nucl_atom_ids[-1][1]+1), range(dna_nucl_sections[-1][0], dna_nucl_sections[-1][1]+1)):
                    bead = dna_full[d_idx]

                    dna_glob.append([i_idx, bead[2]])

                dna_blobs.append(dna_glob)
            
            
            new_blobs = list()
            for glob, dna_glob in zip(C_blobs, dna_blobs):
                new_blobs.append(glob + dna_glob)

            if LH:
                for glob, dna_glob in zip(L_blobs, dna_blobs):
                    new_blobs.append(glob + dna_glob)


            blobs = new_blobs

       
        elif version == 4: # need to add in most of the nucleosomeal dna

            #need to get the atom id of the one before first dna.


            dna_nucl_sections = list()
            dna_nucl_atom_ids = list()
            nrl = int(dna_len_inital / 3)

            mid_point = int(np.floor(nrl / 2))

            start_point = int(mid_point - 147 / 2)
            end_point = int(mid_point + 147 / 2)

            dna_blobs = list()

            print(i)

            for N in range(NUM_NUCLS):
                dna_nucl_sections.append([start_point * 3 + N * nrl * 3, end_point * 3 + N * nrl * 3 + 2])

                dna_nucl_atom_ids.append([dna_nucl_sections[-1][0] + i, dna_nucl_sections[-1][1] + i])

                dna_glob = list()
                for i_idx,d_idx in zip(range(dna_nucl_atom_ids[-1][0], dna_nucl_atom_ids[-1][1]+1), range(dna_nucl_sections[-1][0], dna_nucl_sections[-1][1]+1)):
                    bead = dna_full[d_idx]

                    dna_glob.append([i_idx, bead[2]])

                dna_blobs.append(dna_glob)




            new_blobs = list()
            for glob, dna_glob in zip(C_blobs, dna_blobs):
                new_blobs.append(glob + dna_glob)

            if LH:
                for glob, dna_glob in zip(L_blobs, dna_blobs):
                    new_blobs.append(glob + dna_glob)


            blobs = new_blobs

        else:
            blobs = L_blobs + C_blobs

        #for blob in blobs:
        #    for line in blob:
        #        #print(line)




        # find neighbours of each bead

        bonds_temp = list()
        bond_dist_dict = dict()

        for glob in blobs:
            for bead_i in glob:
                ii = bead_i[0]

                dists = list()
                for bead_j in glob:
                    jj = bead_j[0]
                    if jj != ii and not ( jj >= i  and ii >= i ): # remove and dna-dna bonds
                        dist = np.linalg.norm([bead_i[1] - bead_j[1]])

                        dists.append([jj, dist])

                # find the N smallest
                dists = np.array(dists)
                sorted = dists[dists[:, 1].argsort()]

                for n in range(len(sorted)):
                    if sorted[n][1] < EN_CUT:
                        bi = ii
                        bj = int(sorted[n][0])

                        if bj > bi:

                            bonds_temp.append([bi, bj])

                            bond_dist_dict[str(bi) + "-" + str(bj)] = sorted[n][1]
                        else:
                            bonds_temp.append([bj, bi])
                            bond_dist_dict[str(bj) + "-" + str(bi)] = sorted[n][1]
                    else:
                        break

        # remove duplicate bonds
        new_bonds = list(set(list(map(tuple, bonds_temp))))
        new_bonds.sort()


        t = 4

        for bond in new_bonds:
            dist = bond_dist_dict[str(bond[0]) + "-" + str(bond[1])]
            bonds.append([b, t, bond[0], bond[1]])

            print("bond_coeff " + str(t) + " harmonic 10.0 " + str(dist), file=in_file)

            t = t + 1
            b = b + 1

        if version==5 and LH:
            for bond,dist in zip(dummy_bond_ids,dummy_bond_dist):
                bonds.append([b,t,bond[0],bond[1]])
                print("bond_coeff " + str(t) + " harmonic 0.1 " + str(dist), file=in_file)
                b = b + 1
                t = t + 1




    else:
        for alt_proteins in nucls:
            for protein in alt_proteins:

                s = 1

                for subsection in protein:
                    if s == 1:
                        for atom in subsection:
                            bonds.append([b, 1, i, i + 1])
                            b = b + 1
                            i = i + 1
                    elif s == 2:
                        for atom in subsection:
                            # bonds.append([b, 1, i, i + 1])
                            # b = b + 1
                            i = i + 1
                    elif s == 3:
                        for atom in subsection:
                            bonds.append([b, 1, i - 1, i])
                            b = b + 1
                            i = i + 1

                    s = s + 1

    # bonds for dna, and angles
    angles = list()
    a = 1
    for n in range(1, n_dna):
        bonds.append([b, 2, i, i + 3])  # dna-dna
        b = b + 1

        # make the placeholder bonds for all nieghbouring phosphates so the lammps
        # special fene skips them from the pairwise interactions
        #bonds.append([b,3,i+1,i+2])
        b = b + 1
        #if (n < (n_dna - 1)):
        bonds.append([b,3,i+1,i+4])
        b = b + 1
        bonds.append([b, 3, i + 1, i + 5])
        b = b + 1
        bonds.append([b, 3, i + 2, i + 4])
        b = b + 1
        bonds.append([b, 3, i + 2, i + 5])
        b = b + 1




        # angles.append([a, 1, i, i + 1, i + 3])  # P1
        # a = a + 1
        # angles.append([a, 2, i, i + 2, i + 3])  # P2
        # a = a + 1
        i = i + 3





    #find com and recenter
    R = np.array([0, 0, 0])
    i = 0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
                    R = R + r
                    i = i + 1

    if version == 5 and LH:
        for atom in dummy_atom:
            r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
            R = R + r
            i = i + 1

    for atom in dna_full:
        r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
        R = R + r
        i = i + 1

    R = R / i
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    atom[2][0] = atom[2][0] - R[0]
                    atom[2][1] = atom[2][1] - R[1]
                    atom[2][2] = atom[2][2] - R[2]

    if version == 5 and LH:
        for atom in dummy_atom:
            atom[2][0] = atom[2][0] - R[0]
            atom[2][1] = atom[2][1] - R[1]
            atom[2][2] = atom[2][2] - R[2]

    for atom in dna_full:
        atom[2][0] = atom[2][0] - R[0]
        atom[2][1] = atom[2][1] - R[1]
        atom[2][2] = atom[2][2] - R[2]

    # find box corrds
    max_r = 0.0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
                    max_el = np.max(np.fabs(r))
                    max_r = np.max([max_el, max_r])

    for atom in dna_full:
        r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
        max_el = np.max(np.fabs(r))
        max_r = np.max([max_el, max_r])

    max_r = 1.5 * max_r

    # create lammps input file

    # count total number of beads
    n_beads = 0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    n_beads = n_beads + 1

    if version == 5 and LH:
        for atom in dummy_atom:
            n_beads = n_beads + 1

    n_beads = n_beads + len(dna_full)
    n_bonds = len(bonds)
    n_angles = len(angles)

    print("", file=in_file)
    # print("angle_style zero", file=in_file)
    # print("angle_coeff 1", file=in_file)
    # print("angle_coeff 2", file=in_file)
    # print("", file=in_file)

    # need to work out the pairs coeffs for  the amino acids
    aa_pair_coeff_list = list()

    # we have 40 x 40 combinations

    if HPS:

        for ii in range(1,21):
            for jj in range(ii,21):

                aa_i = index_dict[ii]
                aa_j = index_dict[jj]
                sigma = 0.5*(vdw_dict[aa_i] + vdw_dict[aa_j])

                llambda = 0.5*(hps_lambdas[aa_i] + hps_lambdas[aa_j])

                epsilon = 0.2


                q_1 = charge_dict[aa_i]
                q_2 = charge_dict[aa_j]

                qq = q_1*q_2

                if(qq > 0.0 or qq < 0.0):
                    q_cut = Q_CUT

                else:
                    q_cut = 0.0

                print("pair_coeff " + str(ii) + " " + str(jj) + " " + str(epsilon)
                      + " " + str(sigma) + " " + str(llambda) + " " + str(3*sigma)
                      + " " + str(q_cut), file=in_file)


    else:

        for ii in range(1,41):
            for jj in range(ii,41):

                aa_i = index_dict[ii]
                aa_j = index_dict[jj]
                sigma = 0.5*(vdw_dict[aa_i] + vdw_dict[aa_j])

                KH_key = aa_i + " " + aa_j
                KH_key_rev = aa_j + " " + aa_i

                # if both in tails we use D params:
                if(ii < 21 and jj < 21):



                    try:
                        e = KH_D_dict[KH_key]

                    except KeyError:
                        e = KH_D_dict[KH_key_rev]

                    #epsilon = np.abs(KH_A_a*(e - KH_A_e))

                    # TODO check this, plot graph maybe?
                    #print("ekh = ", e," e0 = ", KH_A_e )
                    # if e < KH_A_e:
                    #     llambda = 1
                    # else:
                    #     llambda = -1

                else:
                    try:
                        e = KH_A_dict[KH_key]

                    except KeyError:
                        e = KH_A_dict[KH_key_rev]

                    #epsilon = np.abs(KH_D_a*(e - KH_D_e))

                    # if e < KH_D_e:
                    #     llambda = 1
                    # else:
                    #     llambda = -1

                epsilon = np.abs(e)
                if e < 0.0:
                    llambda = 1.0
                else:
                    llambda = -1.0


                # if aa_i == "MET" and aa_j == "GLY" or aa_j == "MET" and aa_i == "GLY":
                #     print(epsilon, llambda)


                q_1 = charge_dict[aa_i]
                q_2 = charge_dict[aa_j]

                qq = q_1*q_2

                if(qq > 0.0 or qq < 0.0):
                    q_cut = Q_CUT

                else:
                    q_cut = 0.0


                xx = list(np.linspace(2.0,20.0,1000))

                yy = list()
                for x in xx:
                    yy.append(lk_KA_test_potential(x,sigma,epsilon,llambda))


                #plt.plot(xx,yy)
                #plt.ylim([-1, 1])
                #plt.show()



                print("pair_coeff " + str(ii) + " " + str(jj) + " " + str(epsilon)
                      + " " + str(sigma) + " " + str(llambda) + " " + str(3*sigma)
                      + " " + str(q_cut), file=in_file)

    # compute all pair coeffs for protein dna and protein P
    #plt.show()
    if HPS:

        for key in charge_dict:
            aa = key
            index1 = index_dict_rev[aa]


            aa_vdw = vdw_dict[aa]

            #D_combined_vdw = 0.5*(aa_vdw + DNA_vdw)
            #P_combined_vdw = 0.5*(aa_vdw + P_vdw)

            D_combined_vdw = DNA_vdw
            P_combined_vdw = P_vdw
            str_dna_1 = "pair_coeff " + str(index1)  + " 21 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(D_lambda) + " " + str(3*D_combined_vdw) + " 0.0"
            #str_dna_2 = "pair_coeff " + str(index2)  + " 21 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(D_lambda) + " " + str(3*D_combined_vdw) + " 0.0"

            if charge_dict[aa] > 0.0 or charge_dict[key] < 0.0:
                q_cut = Q_CUT
            else:
                q_cut = 0.0


            str_P_1 = "pair_coeff " + str(index1)  + " 22 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(P_lambda) + " " + str(3*P_combined_vdw) + " " + str(q_cut)
            #str_P_2 = "pair_coeff " + str(index2)  + " 22 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(P_lambda) + " " + str(3*P_combined_vdw) + " " + str(q_cut)

            print(str_dna_1, file=in_file)
            #print(str_dna_2, file=in_file)
            print(str_P_1, file=in_file)
            #print(str_P_2, file=in_file)

        print("pair_coeff         21     21        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         21     22        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         22     22        0.0   0.0      0.0       0.0 " + str(Q_CUT), file=in_file)


    else:
        for key in charge_dict:
            aa = key
            index1 = index_dict_rev[aa]
            index2 = index1 + 20

            aa_vdw = vdw_dict[aa]

            #D_combined_vdw = 0.5*(aa_vdw + DNA_vdw)
            #P_combined_vdw = 0.5*(aa_vdw + P_vdw)

            D_combined_vdw = DNA_vdw
            P_combined_vdw = P_vdw

            str_dna_1 = "pair_coeff " + str(index1)  + " 41 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(D_lambda) + " " + str(3*D_combined_vdw) + " 0.0"
            str_dna_2 = "pair_coeff " + str(index2)  + " 41 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(D_lambda) + " " + str(3*D_combined_vdw) + " 0.0"

            if charge_dict[aa] > 0.0 or charge_dict[key] < 0.0:
                q_cut = Q_CUT
            else:
                q_cut = 0.0


            str_P_1 = "pair_coeff " + str(index1)  + " 42 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(P_lambda) + " " + str(3*P_combined_vdw) + " " + str(q_cut)
            str_P_2 = "pair_coeff " + str(index2)  + " 42 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(P_lambda) + " " + str(3*P_combined_vdw) + " " + str(q_cut)

            print(str_dna_1, file=in_file)
            print(str_dna_2, file=in_file)
            print(str_P_1, file=in_file)
            print(str_P_2, file=in_file)


        print("pair_coeff         41     41        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         41     42        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         42     42        0.0   0.0      0.0       0.0 " + str(Q_CUT), file=in_file)

    if version == 5:
        # extra dummy bead
        if HPS:
            for i in range(1, 21):
                print("pair_coeff   " + str(i) + "           23        0.0   0.0      0.0       0.0     0.000",
                      file=in_file)

            print("pair_coeff   23    23        0.0   0.0      0.0       0.0     00.000", file=in_file)
        else:
            for i in range(1, 41):
                print("pair_coeff   " + str(i) + "           43        0.0   0.0      0.0       0.0     0.000",
                      file=in_file)

            print("pair_coeff   43    43        0.0   0.0      0.0       0.0     00.000", file=in_file)

    # exclude only direct bonded pair interactions
    print("special_bonds fene", file=in_file)
    print("neighbor 10 bin", file=in_file)
    #print("neigh_modify  every 1 delay 0", file=in_file)
    

    if version == 2:
        data_file = open("nucl_v2.txt", "w")
    else:
        data_file = open("nucl.txt", "w")
    
    data_file.write("#LAMMPS data file\n")
    data_file.write("\n")
    data_file.write(str(n_beads) + " atoms\n")
    data_file.write(str(n_dna) + " ellipsoids\n")
    data_file.write(str(n_bonds) + " bonds\n")
    data_file.write(str(n_angles) + " angles\n")
    data_file.write("0 dihedrals\n")
    data_file.write("0 impropers\n")
    data_file.write("\n")

    if version == 5:
        if HPS:
            data_file.write("23 atom types\n")
        else:
            data_file.write("43 atom types\n")

        if ENM:
            data_file.write(str(t - 1) + " bond types\n")
        else:
            data_file.write("3 bond types\n")

    else:

        if HPS:
            data_file.write("22 atom types\n")
        else:
            data_file.write("42 atom types\n")

        if ENM:
            data_file.write(str(t - 1) + " bond types\n")
        else:
            data_file.write("3 bond types\n")
    data_file.write("0 angle types\n")
    data_file.write("\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " xlo xhi\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " ylo yhi\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " zlo zhi\n")

    bp_file = open("DNA_sequence.txt", "w")
    bp_file.write("# " + str(n_dna) + "\n")
    data_file.write("\n")
    data_file.write("Atoms\n")
    data_file.write("\n")
    i = 1
    m = 0

    if version == 1 or version == 3 or version == 4 or version == 5:

        mol_starts = list()
        mol_ends = list()
        for alt_proteins in nucls:
            p = 1
            for protein in alt_proteins:
                if p == 1:
                    m = m + 1
                else:
                    m = m

                mol = m
                s = 1
                for subsection in protein:
                    start_i = i

                    for atom in subsection:
                        a_type = index_dict_rev[atom[0]]

                        if s == 2:
                            # rigid region
                            if not HPS:
                                a_type = a_type + 20

                        charge = charge_dict[atom[0]]
                        # atom id, atom type, x, y, z ellipsoid flag, density,
                        # molecule, charge
                        data_file.write(
                            str(i) + " " + str(a_type) + " " + str(atom[2][0]) + " " +
                            str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                            str(mol) + " " + str(charge) + "\n")

                        i = i + 1
                    end_i = i - 1

                    if s == 2:
                        print(str(start_i) + ":" + str(end_i))
                        mol_starts.append(start_i)
                        mol_ends.append(end_i)
                    s = s + 1
                if LH:
                    if p == 1:
                        m = m + 1
                p = p + 1

            print("\n")
        m = m
        mol = m
        ellipsoid_list = list()
        first_dna_mol = mol+1
        for bead in dna_full:
            # atom id, atom type, x, y, z ellipsoid flag, density, molecule, charge
            if bead[0] == "D":
                a_type = D_type
                el = 1
                charge = 0.0
                ellipsoid_list.append(i)
                bp_file.write(str(i) + " " + bead[1][0] + bead[1][2] + "\n")
                mol = mol + 1
            else:
                a_type = P_type
                el = 0
                charge = -1.0

            data_file.write(
                str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                str(mol) + " " + str(charge) + "\n")

            i = i + 1
        last_dna_mol = mol
    else: # need to group entire nucleosome protein + dna into one molecule
        m=1
        mol_starts = list()
        mol_ends = list()

        for alt_proteins in nucls:
            p = 1
            for protein in alt_proteins:
                # if p == 1:
                #     m = m + 1
                # else:
                #     m = m

                mol = m
                s = 1
                for subsection in protein:
                    start_i = i

                    for atom in subsection:
                        a_type = index_dict_rev[atom[0]]

                        if s == 2:
                            # rigid region
                            if not HPS:
                                a_type = a_type + 20

                        charge = charge_dict[atom[0]]
                        # atom id, atom type, x, y, z ellipsoid flag, density,
                        # molecule, charge
                        data_file.write(
                            str(i) + " " + str(a_type) + " " + str(atom[2][0]) + " " +
                            str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                            str(mol) + " " + str(charge) + "\n")

                        i = i + 1
                    end_i = i - 1

                    if s == 2:
                        print(str(start_i) + ":" + str(end_i))
                        mol_starts.append(start_i)
                        mol_ends.append(end_i)
                    s = s + 1
                # if p == 1:
                #     m = m + 1
                p = p + 1

            m = m + 1

            print("\n")
        m = m
        mol = m

        if version == 5 and LH:
            mmm=1
            for atom in dummy_atom:
                if HPS:
                    X_TYPE=23
                else:
                    X_TYPE=43
                data_file.write(
                            str(i) + " " + str(X_TYPE) + " " + str(atom[2][0]) + " " +
                            str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                            str(mmm) + " " + str(0.0) + "\n")
                i=i+1
                mmm=mmm+1


        ellipsoid_list = list()

        dna_nucl_sections = list()
        dna_nucl_atom_ids = list()
        nrl = int(dna_len_inital/3)

        mid_point = nrl/2

        start_point = int(mid_point-147/2)
        end_point = int(mid_point+147/2)

        for N in range(NUM_NUCLS):
            dna_nucl_sections.append([start_point*3 + N*nrl*3, end_point*3 + N*nrl*3 + 2])

            dna_nucl_atom_ids.append([dna_nucl_sections[-1][0]+i, dna_nucl_sections[-1][1]+i])
        j = 0
        for bead in dna_full:
            # atom id, atom type, x, y, z ellipsoid flag, density, molecule, charge
            if bead[0] == "D":
                a_type = D_type
                el = 1
                charge = 0.0
                ellipsoid_list.append(i)
                bp_file.write(str(i) + " " + bead[1][0] + bead[1][2] + "\n")
                mol = mol + 1
            else:
                a_type = P_type
                el = 0
                charge = -1.0
            n = 1

            written = False
            for nucl_dna in dna_nucl_sections:
                if (j <= nucl_dna[1] and j >= nucl_dna[0]):
                    m = n

                    data_file.write(
                        str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                        str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                        str(m) + " " + str(charge) + "\n")

                    written = True

                n = n + 1

            if not written:
                data_file.write(
                    str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                    str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                    str(mol) + " " + str(charge) + "\n")

            i = i + 1
            j = j + 1

    bp_file.close()

    # if seq_type == 5:
    #     seq_file = open("DNA_sequence.txt",'r')
    #     data = list()
    #     for line in seq_file:
    #         data.append(line)
    #
    #     seq_file.close()
    #     s_len = len(data)
    #
    #     data_array = list()
    #     for k in range(1,s_len):
    #         split_s = data[k].split()
    #
    #         data_array.append([split_s[0], split_s[1]])
    #
    #     data_array = np.array(data_array)
    #     diff = new_nrl_length - 211
    #     ncp147 = data_array[int(np.floor(diff/2)):-int(np.ceil(diff/2)),1]
    #
    #
    #
    #     v1 = data_array[:]
    #     v1[:211,1] = ncp147[:]
    #     v1[-211:,1] = ncp147[:]
    #     v1[211:-211] = "AT"
    #
    #     v2 = v1[:]
    #
    #     v2[211:-211] = "TA"
    #
    #     v1_file = open("DNA_sequence_2.txt",'w')
    #     v2_file = open("DNA_sequence_3.txt", 'w')
    #     v1_file.write(data[0])
    #     v2_file.write(data[0])
    #
    #     for k in range(0,s_len-1):
    #         v1_file.write(v1[k,0] + " " + v1[k,1]+"\n")
    #         v2_file.write(v2[k,0] + " " + v2[k,1]+"\n")
    #     v1_file.close()
    #     v2_file.close()

    if seq_type == 5:
        seq_file = open("DNA_sequence.txt",'r')
        data = list()
        for line in seq_file:
            data.append(line)

        seq_file.close()
        s_len = len(data)

        col_1=list()
        col_2=list()


        for k in range(1,s_len):
            split_s = data[k].split()

            col_1.append(split_s[0])
            col_2.append(split_s[1])



        v1_file = open("DNA_sequence_2.txt",'w')
        v1_file.write(data[0])

        new_col_2 = list()
        for line in col_2:
            new_col_2.append(line)

        for i in range(left_l,left_l+211):
            new_col_2[i] = col_2[-(211+left_l)+(i-left_l)]
            new_col_2[-(211+left_l)+(i-left_l)] = col_2[i]

        for k in range(0,s_len-1):
            v1_file.write(col_1[k] + " " + new_col_2[k]+"\n")
        v1_file.close()

    nucl_string = "group nucl id"
    for start, end in zip(mol_starts, mol_ends):
        nucl_string = nucl_string + " " + str(start) + ":" + str(end)

    if version == 2:
        print("here")
        for dna_bloc in dna_nucl_atom_ids:
            nucl_string = nucl_string + " " + str(dna_bloc[0]) + ":" + str(dna_bloc[1])

    print(nucl_string, file=in_file)

    data_file.write("\n")
    data_file.write("Ellipsoids\n")
    data_file.write("\n")
    #print("group dna_all id " + str(ellipsoid_list[0]) + ":" + str(ellipsoid_list[-1]), file=in_file)
    #print("group notdna subtract all dna_all", file=in_file)
    # print("group tails  subtract notdna nucl", file=in_file)
    #print("group dna id " + str(ellipsoid_list[0]) + ":" + str(ellipsoid_list[-1]) + ":3", file=in_file)
    a = 0
    for bead in dna_full:
        if bead[0] == "D":
            data_file.write(
                str(ellipsoid_list[a]) + " " + str(RX) + " " + str(RY) + " " +
                str(RZ) + " " + str(bead[3][0]) + " " + str(bead[3][1]) + " " +
                str(bead[3][2]) + " " + str(bead[3][3]) + "\n")
            a = a + 1

    data_file.write("\n")
    data_file.write("Bonds\n")
    data_file.write("\n")
    for bond in bonds:
        data_file.write(
            str(bond[0]) + " " + str(bond[1]) + " " + str(bond[2]) + " " +
            str(bond[3]) + "\n")

    data_file.write("\n")
    # data_file.write("Angles\n")
    # data_file.write("\n")
    # for angle in angles:
    #     data_file.write(
    #         str(angle[0]) + " " + str(angle[1]) + " " + str(angle[2]) + " " +
    #         str(angle[3]) + " " + str(angle[4]) + "\n")
    # data_file.close()

    # print("group end1 id " + str(ellipsoid_list[0]), file=in_file)
    # print("group end2 id " + str(ellipsoid_list[-1]), file=in_file)
    #print("group dna_P subtract dna_all dna", file=in_file)

    if HPS:
        print("group dna type 21", file=in_file)
        print("group ps type 22", file=in_file)
    else:
        print("group dna type 41", file=in_file)
        print("group ps type 42", file=in_file)


    if version == 1 or version == 3 or version == 4 or version == 5:
        print("group dna_all union dna ps", file=in_file)
        print("group notdna subtract all dna_all", file=in_file)

        print("neigh_modify exclude molecule/intra nucl", file=in_file)
        print("neigh_modify exclude molecule/intra dna_all ", file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(P_type), file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(D_type), file=in_file)

        first_dna_id = ellipsoid_list[0]
        last_dna_id = ellipsoid_list[-1]
        
        #nrl = (last_dna_id - first_dna_id)/3
        #diff = nrl - 211
        #half = int(diff/2)*3
        half =0 
        print('variable ex equal "sqrt((x['+str(first_dna_id+half)+'] - x['+str(last_dna_id-half)+'])^2 + (y['+str(first_dna_id+half)+'] - y['+str(last_dna_id-half)+'])^2 + (z['+str(first_dna_id+half)+'] - z['+str(last_dna_id-half)+'])^2)"',file=in_file)
        #print('group end1 molecule '+str(first_dna_mol),file=in_file)
        #print('group end2 molecule '+str(last_dna_mol),file=in_file)
        print('group end1 id '+str(first_dna_id),file=in_file)
        print('group end2 id '+str(last_dna_id),file=in_file)
        print('group P1 id ' + str(first_dna_id + 1),file=in_file)
        print('group P2 id ' + str(first_dna_id + 2),file=in_file)
        print('group p3 id ' + str(last_dna_id  + 1),file=in_file)
        print('group p4 id ' + str(last_dna_id  + 2),file=in_file)
        
        # if version == 4:
        #     settings = open("lmp_run_settings_nucl_bead_string_setup.txt", 'r')
        #     for line in settings:
        #         in_file.write(line)
        # else:
        #     if not RESTART:
        #         settings = open("lmp_run_settings_nucl.txt", 'r')
        #         for line in settings:
        #             if not version == 5:
        #                 if line == "comm_modify cutoff 80\n":
        #                     line="\n"

        #             in_file.write(line)
        #     else:
        #         settings = open("lmp_run_settings_nucl_restart.txt", 'r')
        #         for line in settings:
        #             if not version == 5:
        #                 if line == "comm_modify cutoff 80\n":
        #                     line = "\n"
        #             in_file.write(line)
    else:
        print("group dna_all_temp union dna ps", file=in_file)
        print("group dna_all subtract dna_all_temp nucl", file=in_file)
        print("group temp subtract all nucl", file=in_file)
        print("group tails subtract temp dna_all", file=in_file)
        print("group dna_only intersect dna_all dna", file=in_file)


        print("neigh_modify exclude molecule/intra nucl", file=in_file)
        print("neigh_modify exclude molecule/intra dna_all ", file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(P_type), file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(D_type), file=in_file)


        settings = open("lmp_run_settings_v2.txt", 'r')
        for line in settings:
            in_file.write(line)




    in_file.close()


def create_nucl_from_dump(EPSILON_P,EPSILON_D, HPS,NUM_NUCLS, version, LH, seq_type,RESTART,ADD_NRL,left_l, right_l,RESTN,ishremd):
    # ------------------------------------------------------------------------------
    # Dictionaries for amino acid parameters

    # 601 sequence:
    seq601 = "ACAGGATGTATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGACAGCGCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTCTCCAG"
    seqTTAGGG = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA"
    seqAT = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    RISE = 3.4
    TWIST = 36 * np.pi / 180.0

    # Charges
    charge_dict = {
        "ALA": 0,
        "ARG": 1,
        "ASN": 0,
        "ASP": -1,
        "CYS": 0,
        "GLN": 0,
        "GLU": -1,
        "GLY": 0,
        "HIS": 0.5,
        "ILE": 0,
        "LEU": 0,
        "LYS": 1,
        "MET": 0,
        "PHE": 0,
        "PRO": 0,
        "SER": 0,
        "THR": 0,
        "TRP": 0,
        "TYR": 0,
        "VAL": 0,
    }

    # VdW radii
    vdw_dict = {
        "ALA": 5.04,
        "ARG": 6.56,
        "ASN": 5.68,
        "ASP": 5.58,
        "CYS": 5.48,
        "CYD": 5.48,
        "GLN": 6.02,
        "GLU": 5.92,
        "GLY": 4.50,
        "HIS": 6.08,
        "ILE": 6.18,
        "LEU": 6.18,
        "LYS": 6.36,
        "MET": 6.18,
        "PHE": 6.36,
        "PRO": 5.56,
        "SER": 5.18,
        "THR": 5.62,
        "TRP": 6.78,
        "TYR": 6.46,
        "VAL": 5.86,
    }

    hps_lambdas = {
        "ILE": 0.972973,
        "VAL": 0.891892,
        "LEU": 0.972973,
        "PHE": 1.0,
        "CYS": 0.594595,
        "MET": 0.837838,
        "ALA": 0.72973,
        "GLY": 0.648649,
        "THR": 0.675676,
        "SER": 0.594595,
        "TRP": 0.945946,
        "TYR": 0.864865,
        "PRO": 1.0,
        "HIS": 0.513514,
        "GLU": 0.459459,
        "GLN": 0.513514,
        "ASP": 0.378378,
        "ASN": 0.432432,
        "LYS": 0.513514,
        "ARG": 0.0
    }

    # masses
    mass_dict = {
        "ALA": 71.08,
        "ARG": 156.20,
        "ASN": 114.10,
        "ASP": 115.10,
        "CYS": 103.10,
        "GLN": 128.10,
        "GLU": 129.10,
        "GLY": 57.05,
        "HIS": 137.10,
        "ILE": 113.20,
        "LEU": 113.20,
        "LYS": 128.20,
        "MET": 131.20,
        "PHE": 147.20,
        "PRO": 97.12,
        "SER": 87.08,
        "THR": 101.10,
        "TRP": 186.20,
        "TYR": 163.20,
        "VAL": 99.07,
    }

    """ Main program --------------------------------------------------------"""

    D_a_epsilon = EPSILON_D
    P_a_epsilon = EPSILON_P

    P_lambda = 1.0
    D_lambda = 1.0

    # use elastic network model

    if version == 1 or version == 3 or version == 4 or version == 5:

        ENM = True
    else:
        ENM = False

    # cut-off distance for creating the elastic network
    EN_CUT = 7.5
    HPS_MOD = False
    if HPS_MOD:
        HPS = False

    if HPS:
        D_type = 21
        P_type = 22

    else:
        D_type = 41
        P_type = 42

    DNA_vdw = 8
    P_vdw = 4

    D_a_sigma_dict = {}

    P_a_sigma_dict = {}

    WIDTH = 8.5
    RX = 11.0
    RY = 20.0
    RZ = 3.5

    P_ANGLE = 20.0 * np.pi / 180.0

    #Q_CUT = 3.5*DEBYE_LENGTH
    Q_CUT = "$Q"

    if not HPS:
        # get the KH params dicts
        file_A = open("KH_A.txt", "r")

        KH_A_dict = {}

        for line in file_A:
            split = line.split()
            KH_A_dict[split[0] + " " + split[1]] = float(split[2])

        file_D = open("KH_D.txt", "r")

        KH_D_dict = {}

        for line in file_D:
            split = line.split()
            KH_D_dict[split[0] + " " + split[1]] = float(split[2])

    if version == 2:
        in_file = open("in.nucl_v2", "w")
    else:
        if not RESTART:
            in_file = open("in.nucl", "w")
        else:
            in_file = open("in.nucl_restart", "w")

    # VARIABLES
    #if version == 2:
    #    print("variable fname index nucl_v2.txt", file=in_file)
    #    print("variable simname index nucl_v2", file=in_file)
    #else:
    #    print("variable fname index nucl_v1.txt", file=in_file)
    #    print("variable simname index nucl_v1", file=in_file)

    print("", file=in_file)
    print("", file=in_file)
    print("variable a equal 8.0", file=in_file)
    print("variable l equal 1.0/${a}", file=in_file)
    print("variable Q equal 3.5*${a}", file=in_file)
    print("", file=in_file)
    #print("variable I world A B", file=in_file)
    print("# Initialization", file=in_file)
    print("units		real", file=in_file)
    print("atom_style   hybrid ellipsoid angle charge", file=in_file)
    print("", file=in_file)
    print("boundary s s s", file=in_file)
    print("", file=in_file)
    if not RESTART:
        if ishremd:
            print("log log.temper.${I}.txt", file=in_file)
        else:
            print("log log.nucl.txt", file=in_file)

        
        print("read_data nucl.txt", file=in_file)
    else:
        print("log 	log.${simname}.txt append", file=in_file)
        print("read_restart nucl.restart.*", file=in_file)

    print("", file=in_file)

    aa_keys = list(mass_dict)
    aa_keys.sort()

    # Masses
    # "set type    X mass  M"
    index_dict = {}
    index_dict_rev = {}
    a = 1
    for aa in aa_keys:
        print("set type " + str(a) + " mass " + str(mass_dict[aa]), file=in_file)
        index_dict[a] = aa
        index_dict_rev[aa] = a
        a = a + 1

    if not HPS:
        for aa in aa_keys:
            print("set type " + str(a) + " mass " + str(mass_dict[aa]), file=in_file)
            index_dict[a] = aa
            a = a + 1

    print("set type " + str(D_type) + " mass 650.0", file=in_file)
    print("set type " + str(P_type) + " mass 0.000001", file=in_file)

    if version == 5 and LH:
        if HPS:
            print("set type " + str(23) + " mass 100.0", file=in_file)
        else:
            print("set type " + str(43) + " mass 100.0", file=in_file)

    #print("pair_style ljlambda 0.125 " + str(Q_CUT) + " " + str(Q_CUT), file=in_file)
    #print("pair_style ljlambda "+str(1.0/DEBYE_LENGTH) +" " + str(Q_CUT) + " " + str(Q_CUT), file=in_file)
    print("pair_style ljlambda $l 20.34 $Q", file=in_file)
    print("dielectric  80.0", file=in_file)
    print("", file=in_file)
    print("bond_style hybrid harmonic harmonic/DNA zero", file=in_file)
    print("", file=in_file)
    print("bond_coeff 1 harmonic 10.0 3.8", file=in_file)
    print("bond_coeff 2 harmonic/DNA 0.0 0.0", file=in_file)
    print("bond_coeff 3 zero", file=in_file)

    # get_min_distance_mean_all("clust.pdb")
    # get_min_distance_CAs("clust.pdb")

    # get the single nucleosome structure
    # DNA and histones
    alt_proteins, dna_full, n_dna = get_single_nucl(WIDTH, P_ANGLE)

    # load in the dump file and edit to be the new positions
    if HPS:
        D_ID = 21
        P_ID = 22
    else:
        D_ID = 41
        P_ID = 42

    alt_proteins, dna_full, n_dna = update_from_dump(alt_proteins, dna_full, n_dna, "bead_string_setup_last_frame.dump",
                                                     [D_ID, P_ID])

    if ADD_NRL:
        dna_full = add_nrl(dna_full, left_l, right_l, WIDTH, P_ANGLE)

    dna_len_inital = len(dna_full)

    current_nrl = dna_len_inital / 3

    # replace the dna sequence
    if seq_type == 2:
        # only replace inner 147
        start_index = (current_nrl - 147) / 2
        end_index = current_nrl - start_index

        print(start_index)
        print(end_index)

        for iiii in range(147):
            base = seq601[iiii]
            if base == "T":
                base2 = "A"
            elif base == "A":
                base2 = "T"
            elif base == "C":
                base2 = "G"
            else:
                base2 = "C"

            dna_idx = int((iiii + start_index) * 3)
            dna_full[dna_idx][1] = base + "-" + base2

    if seq_type == 3:
        # only replace inner 147
        start_index = (current_nrl - 147) / 2
        end_index = current_nrl - start_index

        print(start_index)
        print(end_index)

        for iiii in range(147):
            base = seqTTAGGG[iiii]
            if base == "T":
                base2 = "A"
            elif base == "A":
                base2 = "T"
            elif base == "C":
                base2 = "G"
            else:
                base2 = "C"

            dna_idx = int((iiii + start_index) * 3)
            dna_full[dna_idx][1] = base + "-" + base2

    if seq_type == 4:
        # replace all
        for iiii in range(int(len(dna_full) / 3)):
            dna_idx = int(iiii * 3)
            dna_full[dna_idx][1] = "A-T"

    if seq_type == "random":
        for iiii in range(int(len(dna_full) / 3)):
            dna_idx = int(iiii * 3)
            dna_full[dna_idx][1] = np.random.choice(["A-T", "T-A", "C-G", "G-C"])

    if not LH:
        alt_proteins = alt_proteins[1:]

    nucls = [alt_proteins]

    for N in range(1, NUM_NUCLS):
        # now we add in an extra nucleosomes
        new_alt_proteins, dna_full = add_nucl(alt_proteins, dna_full, dna_len_inital)

        nucls = nucls + [new_alt_proteins]

    # count number of DNA basepairs
    n_dna = 0
    for dna in dna_full:
        if dna[0] == "D":
            n_dna = n_dna + 1

    # sort out the bonds
    bonds = list()

    # counters for bond and bead index
    b = 1
    i = 1

    if ENM:
        # more complex process required here
        # need to make a network of all globular beads
        # need to identify their nearest neighbours

        # first assign tails as normal, and group the globs into seperate structure

        C_blobs = list()
        L_blobs = list()
        for alt_proteins in nucls:
            core_glob = list()
            H1_glob = list()
            p = 1
            for protein in alt_proteins:

                s = 1

                for subsection in protein:
                    if s == 1:
                        for atom in subsection:
                            bonds.append([b, 1, i, i + 1])
                            b = b + 1
                            i = i + 1
                    elif s == 2:
                        for atom in subsection:
                            if LH:
                                if p == 1:
                                    H1_glob.append([i, atom[2]])
                                else:
                                    core_glob.append([i, atom[2]])
                            else:
                                core_glob.append([i, atom[2]])

                            i = i + 1


                    elif s == 3:
                        for atom in subsection:
                            bonds.append([b, 1, i - 1, i])
                            b = b + 1
                            i = i + 1

                    s = s + 1

                p = p + 1
            C_blobs.append(core_glob)
            L_blobs.append(H1_glob)

        if version == 5 and LH:

            # get the com and add in dummy bead
            C_coms = list()
            L_coms = list()

            for blob in L_blobs:
                pos = list()
                for atom in blob:
                    pos.append(atom[1])

                pos = np.array(pos)
                com = np.mean(pos, axis=0)
                #print(com)
                L_coms.append(com)
            for blob in C_blobs:
                pos = list()
                for atom in blob:
                    pos.append(atom[1])

                pos = np.array(pos)
                com = np.mean(pos, axis=0)
                #print(com)
                C_coms.append(com)
            dummy_bond_ids = list()
            dummy_bond_dist = list()
            dummy_atom = list()
            for cc in range(len(C_blobs)):
                L_blobs[cc].append([i, L_coms[cc]])
                dummy_atom.append(["X", "X", L_coms[cc]])
                i = i + 1
                C_blobs[cc].append([i, C_coms[cc]])
                dummy_atom.append(["X", "X", C_coms[cc]])
                dummy_bond_ids.append([i, i - 1])
                dummy_bond_dist.append(np.linalg.norm(C_coms[cc] - L_coms[cc]))
                i = i + 1

        # print(core_glob)

        if version == 3:  # need to add in most of the nucleosomeal dna

            # need to get the atom id of the one before first dna.

            dna_nucl_sections = list()
            dna_nucl_atom_ids = list()
            nrl = int(dna_len_inital / 3)

            mid_point = int(np.floor(nrl / 2))
            
            #TODO: amount unwrapped
            start_point = int(mid_point - RESTN / 2)
            end_point = int(mid_point + RESTN / 2)

            dna_blobs = list()

            #print(i)

            for N in range(NUM_NUCLS):
                dna_nucl_sections.append([start_point * 3 + N * nrl * 3, end_point * 3 + N * nrl * 3 + 2])

                dna_nucl_atom_ids.append([dna_nucl_sections[-1][0] + i, dna_nucl_sections[-1][1] + i])

                dna_glob = list()
                for i_idx, d_idx in zip(range(dna_nucl_atom_ids[-1][0], dna_nucl_atom_ids[-1][1] + 1),
                                        range(dna_nucl_sections[-1][0], dna_nucl_sections[-1][1] + 1)):
                    bead = dna_full[d_idx]

                    dna_glob.append([i_idx, bead[2]])

                dna_blobs.append(dna_glob)

            new_blobs = list()
            for glob, dna_glob in zip(C_blobs, dna_blobs):
                new_blobs.append(glob + dna_glob)

            if LH:
                for glob, dna_glob in zip(L_blobs, dna_blobs):
                    new_blobs.append(glob + dna_glob)

            blobs = new_blobs


        elif version == 4:  # need to add in most of the nucleosomeal dna

            # need to get the atom id of the one before first dna.

            dna_nucl_sections = list()
            dna_nucl_atom_ids = list()
            nrl = int(dna_len_inital / 3)

            mid_point = int(np.floor(nrl / 2))

            start_point = int(mid_point - 136 / 2)
            end_point = int(mid_point + 136 / 2)

            dna_blobs = list()

            #print(i)

            for N in range(NUM_NUCLS):
                dna_nucl_sections.append([start_point * 3 + N * nrl * 3, end_point * 3 + N * nrl * 3 + 2])

                dna_nucl_atom_ids.append([dna_nucl_sections[-1][0] + i, dna_nucl_sections[-1][1] + i])

                dna_glob = list()
                for i_idx, d_idx in zip(range(dna_nucl_atom_ids[-1][0], dna_nucl_atom_ids[-1][1] + 1),
                                        range(dna_nucl_sections[-1][0], dna_nucl_sections[-1][1] + 1)):
                    bead = dna_full[d_idx]

                    dna_glob.append([i_idx, bead[2]])

                dna_blobs.append(dna_glob)

            new_blobs = list()
            for glob, dna_glob in zip(C_blobs, dna_blobs):
                new_blobs.append(glob + dna_glob)

            if LH:
                for glob, dna_glob in zip(L_blobs, dna_blobs):
                    new_blobs.append(glob + dna_glob)

            blobs = new_blobs

        else:
            blobs = L_blobs + C_blobs

        for blob in blobs:
            for line in blob:
                #print(line)
                pass

        # find neighbours of each bead

        bonds_temp = list()
        bond_dist_dict = dict()

        for glob in blobs:
            for bead_i in glob:
                ii = bead_i[0]

                dists = list()
                for bead_j in glob:
                    jj = bead_j[0]
                    if jj != ii and not (jj >= i and ii >= i):
                        dist = np.linalg.norm([bead_i[1] - bead_j[1]])

                        dists.append([jj, dist])

                # find the N smallest
                dists = np.array(dists)
                sorted = dists[dists[:, 1].argsort()]

                for n in range(len(sorted)):
                    if sorted[n][1] < EN_CUT:
                        bi = ii
                        bj = int(sorted[n][0])

                        if bj > bi:

                            bonds_temp.append([bi, bj])

                            bond_dist_dict[str(bi) + "-" + str(bj)] = sorted[n][1]
                        else:
                            bonds_temp.append([bj, bi])
                            bond_dist_dict[str(bj) + "-" + str(bi)] = sorted[n][1]
                    else:
                        break

        # remove duplicate bonds
        new_bonds = list(set(list(map(tuple, bonds_temp))))
        new_bonds.sort()

        t = 4

        for bond in new_bonds:
            dist = bond_dist_dict[str(bond[0]) + "-" + str(bond[1])]
            bonds.append([b, t, bond[0], bond[1]])

            print("bond_coeff " + str(t) + " harmonic 10.0 " + str(dist), file=in_file)

            t = t + 1
            b = b + 1

        if version == 5 and LH:
            for bond, dist in zip(dummy_bond_ids, dummy_bond_dist):
                bonds.append([b, t, bond[0], bond[1]])
                # igore computed dist, we know what it is
                dist=62
                print("bond_coeff " + str(t) + " harmonic 0.1 " + str(dist), file=in_file)
                b = b + 1
                t = t + 1



    else:
        for alt_proteins in nucls:
            for protein in alt_proteins:

                s = 1

                for subsection in protein:
                    if s == 1:
                        for atom in subsection:
                            bonds.append([b, 1, i, i + 1])
                            b = b + 1
                            i = i + 1
                    elif s == 2:
                        for atom in subsection:
                            # bonds.append([b, 1, i, i + 1])
                            # b = b + 1
                            i = i + 1
                    elif s == 3:
                        for atom in subsection:
                            bonds.append([b, 1, i - 1, i])
                            b = b + 1
                            i = i + 1

                    s = s + 1

    # bonds for dna, and angles
    angles = list()
    a = 1
    for n in range(1, n_dna):
        bonds.append([b, 2, i, i + 3])  # dna-dna
        b = b + 1

        # make the placeholder bonds for all nieghbouring phosphates so the lammps
        # special fene skips them from the pairwise interactions
        # bonds.append([b,3,i+1,i+2])
        b = b + 1
        # if (n < (n_dna - 1)):
        bonds.append([b, 3, i + 1, i + 4])
        b = b + 1
        bonds.append([b, 3, i + 1, i + 5])
        b = b + 1
        bonds.append([b, 3, i + 2, i + 4])
        b = b + 1
        bonds.append([b, 3, i + 2, i + 5])
        b = b + 1

        # angles.append([a, 1, i, i + 1, i + 3])  # P1
        # a = a + 1
        # angles.append([a, 2, i, i + 2, i + 3])  # P2
        # a = a + 1
        i = i + 3

    # find com and recenter
    R = np.array([0, 0, 0])
    i = 0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
                    R = R + r
                    i = i + 1

    if version == 5 and LH:
        for atom in dummy_atom:
            r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
            R = R + r
            i = i + 1

    for atom in dna_full:
        r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
        R = R + r
        i = i + 1

    R = R / i
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    atom[2][0] = atom[2][0] - R[0]
                    atom[2][1] = atom[2][1] - R[1]
                    atom[2][2] = atom[2][2] - R[2]

    if version == 5 and LH:
        for atom in dummy_atom:
            atom[2][0] = atom[2][0] - R[0]
            atom[2][1] = atom[2][1] - R[1]
            atom[2][2] = atom[2][2] - R[2]

    for atom in dna_full:
        atom[2][0] = atom[2][0] - R[0]
        atom[2][1] = atom[2][1] - R[1]
        atom[2][2] = atom[2][2] - R[2]

    # find box corrds
    max_r = 0.0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
                    max_el = np.max(np.fabs(r))
                    max_r = np.max([max_el, max_r])

    for atom in dna_full:
        r = np.array([atom[2][0], atom[2][1], atom[2][2]], dtype=float)
        max_el = np.max(np.fabs(r))
        max_r = np.max([max_el, max_r])

    max_r = 1.5 * max_r

    # create lammps input file

    # count total number of beads
    n_beads = 0
    for alt_proteins in nucls:
        for protein in alt_proteins:
            for subsection in protein:
                for atom in subsection:
                    n_beads = n_beads + 1

    if version == 5 and LH:
        for atom in dummy_atom:
            n_beads = n_beads + 1

    n_beads = n_beads + len(dna_full)
    n_bonds = len(bonds)
    n_angles = len(angles)

    print("", file=in_file)
    # print("angle_style zero", file=in_file)
    # print("angle_coeff 1", file=in_file)
    # print("angle_coeff 2", file=in_file)
    # print("", file=in_file)

    # need to work out the pairs coeffs for  the amino acids
    aa_pair_coeff_list = list()

    # we have 40 x 40 combinations

    if HPS:

        for ii in range(1, 21):
            for jj in range(ii, 21):

                aa_i = index_dict[ii]
                aa_j = index_dict[jj]
                sigma = 0.5 * (vdw_dict[aa_i] + vdw_dict[aa_j])

                llambda = 0.5 * (hps_lambdas[aa_i] + hps_lambdas[aa_j])

                epsilon = 0.2

                q_1 = charge_dict[aa_i]
                q_2 = charge_dict[aa_j]

                qq = q_1 * q_2

                if (qq > 0.0 or qq < 0.0):
                    q_cut = Q_CUT

                else:
                    q_cut = 0.0

                print("pair_coeff " + str(ii) + " " + str(jj) + " " + str(epsilon)
                      + " " + str(sigma) + " " + str(llambda) + " " + str(3 * sigma)
                      + " " + str(q_cut), file=in_file)


    elif HPS_MOD:
        for ii in range(1, 41):
            for jj in range(ii, 41):
                if ii > 20:
                    iii = ii - 20
                else:
                    iii = ii
                if jj > 20:
                    jjj = jj - 20
                else:
                    jjj = jj
                aa_i = index_dict[iii]
                aa_j = index_dict[jjj]
                sigma = 0.5 * (vdw_dict[aa_i] + vdw_dict[aa_j])

                llambda = 0.5 * (hps_lambdas[aa_i] + hps_lambdas[aa_j])
                if ii > 20 and jj > 20:
                    llambda = 0.0

                epsilon = 0.2

                q_1 = charge_dict[aa_i]
                q_2 = charge_dict[aa_j]

                qq = q_1 * q_2

                if (qq > 0.0 or qq < 0.0):
                    q_cut = Q_CUT

                else:
                    q_cut = 0.0

                print("pair_coeff " + str(ii) + " " + str(jj) + " " + str(epsilon)
                      + " " + str(sigma) + " " + str(llambda) + " " + str(3 * sigma)
                      + " " + str(q_cut), file=in_file)



    else:
        for ii in range(1, 41):
            for jj in range(ii, 41):

                aa_i = index_dict[ii]
                aa_j = index_dict[jj]
                sigma = 0.5 * (vdw_dict[aa_i] + vdw_dict[aa_j])

                KH_key = aa_i + " " + aa_j
                KH_key_rev = aa_j + " " + aa_i

                # if both in tails we use D params:
                if (ii < 21 and jj < 21):
                    try:
                        e = KH_D_dict[KH_key]

                    except KeyError:
                        e = KH_D_dict[KH_key_rev]

                    # epsilon = np.abs(KH_A_a*(e - KH_A_e))

                    # TODO check this, plot graph maybe?
                    # print("ekh = ", e," e0 = ", KH_A_e )
                    # if e < KH_A_e:
                    #     llambda = 1
                    # else:
                    #     llambda = -1

                else:
                    try:
                        e = KH_A_dict[KH_key]

                    except KeyError:
                        e = KH_A_dict[KH_key_rev]

                    # epsilon = np.abs(KH_D_a*(e - KH_D_e))

                    # if e < KH_D_e:
                    #     llambda = 1
                    # else:
                    #     llambda = -1

                epsilon = np.abs(e)
                if e < 0.0:
                    llambda = 1.0
                else:
                    llambda = -1.0

                # if aa_i == "MET" and aa_j == "GLY" or aa_j == "MET" and aa_i == "GLY":
                #     print(epsilon, llambda)

                q_1 = charge_dict[aa_i]
                q_2 = charge_dict[aa_j]

                qq = q_1 * q_2

                if (qq > 0.0 or qq < 0.0):
                    q_cut = Q_CUT

                else:
                    q_cut = 0.0

                xx = list(np.linspace(2.0, 20.0, 1000))

                yy = list()
                for x in xx:
                    yy.append(lk_KA_test_potential(x, sigma, epsilon, llambda))

                # plt.plot(xx,yy)
                # plt.ylim([-1, 1])
                # plt.show()

                print("pair_coeff " + str(ii) + " " + str(jj) + " " + str(epsilon)
                      + " " + str(sigma) + " " + str(llambda) + " " + str(3 * sigma)
                      + " " + str(q_cut), file=in_file)

    # compute all pair coeffs for protein dna and protein P
    # plt.show()
    if HPS:

        for key in charge_dict:
            aa = key
            index1 = index_dict_rev[aa]

            aa_vdw = vdw_dict[aa]

            # D_combined_vdw = 0.5*(aa_vdw + DNA_vdw)
            # P_combined_vdw = 0.5*(aa_vdw + P_vdw)

            D_combined_vdw = DNA_vdw
            P_combined_vdw = P_vdw
            str_dna_1 = "pair_coeff " + str(index1) + " 21 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(
                D_lambda) + " " + str(3 * D_combined_vdw) + " 0.0"
            # str_dna_2 = "pair_coeff " + str(index2)  + " 21 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(D_lambda) + " " + str(3*D_combined_vdw) + " 0.0"

            if charge_dict[aa] > 0.0 or charge_dict[key] < 0.0:
                q_cut = Q_CUT
            else:
                q_cut = 0.0

            str_P_1 = "pair_coeff " + str(index1) + " 22 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(
                P_lambda) + " " + str(3 * P_combined_vdw) + " " + str(q_cut)
            # str_P_2 = "pair_coeff " + str(index2)  + " 22 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(P_lambda) + " " + str(3*P_combined_vdw) + " " + str(q_cut)

            print(str_dna_1, file=in_file)
            # print(str_dna_2, file=in_file)
            print(str_P_1, file=in_file)
            # print(str_P_2, file=in_file)

        print("pair_coeff         21     21        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         21     22        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         22     22        0.0   0.0      0.0       0.0 " + str(Q_CUT), file=in_file)


    else:
        for key in charge_dict:
            aa = key
            index1 = index_dict_rev[aa]
            index2 = index1 + 20

            aa_vdw = vdw_dict[aa]

            # D_combined_vdw = 0.5*(aa_vdw + DNA_vdw)
            # P_combined_vdw = 0.5*(aa_vdw + P_vdw)

            D_combined_vdw = DNA_vdw
            P_combined_vdw = P_vdw

            str_dna_1 = "pair_coeff " + str(index1) + " 41 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(
                D_lambda) + " " + str(3 * D_combined_vdw) + " 0.0"
            str_dna_2 = "pair_coeff " + str(index2) + " 41 " + str(D_a_epsilon) + " " + str(D_combined_vdw) + " " + str(
                D_lambda) + " " + str(3 * D_combined_vdw) + " 0.0"

            if charge_dict[aa] > 0.0 or charge_dict[key] < 0.0:
                q_cut = Q_CUT
            else:
                q_cut = 0.0

            str_P_1 = "pair_coeff " + str(index1) + " 42 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(
                P_lambda) + " " + str(3 * P_combined_vdw) + " " + str(q_cut)
            str_P_2 = "pair_coeff " + str(index2) + " 42 " + str(P_a_epsilon) + " " + str(P_combined_vdw) + " " + str(
                P_lambda) + " " + str(3 * P_combined_vdw) + " " + str(q_cut)

            print(str_dna_1, file=in_file)
            print(str_dna_2, file=in_file)
            print(str_P_1, file=in_file)
            print(str_P_2, file=in_file)

        print("pair_coeff         41     41        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         41     42        0.0   0.0      0.0       0.0     0.000", file=in_file)
        print("pair_coeff         42     42        0.0   0.0      0.0       0.0 " + str(Q_CUT), file=in_file)

    if version == 5:
        # extra dummy bead
        if HPS:
            for i in range(1, 21):
                print("pair_coeff   " + str(i) + "           23        0.0   0.0      0.0       0.0     0.000",
                      file=in_file)

            print("pair_coeff   23    23        0.0   0.0      0.0       0.0     00.000", file=in_file)
        else:
            for i in range(1, 41):
                print("pair_coeff   " + str(i) + "           43        0.0   0.0      0.0       0.0     0.000",
                      file=in_file)

            print("pair_coeff   43    43        0.0   0.0      0.0       0.0     00.000", file=in_file)

    # exclude only direct bonded pair interactions
    print("special_bonds fene", file=in_file)
    print("neighbor 10 bin", file=in_file)
    # print("neigh_modify  every 1 delay 0", file=in_file)

    if version == 2:
        data_file = open("nucl_v2.txt", "w")
    else:
        data_file = open("nucl.txt", "w")

    data_file.write("#LAMMPS data file\n")
    data_file.write("\n")
    data_file.write(str(n_beads) + " atoms\n")
    data_file.write(str(n_dna) + " ellipsoids\n")
    data_file.write(str(n_bonds) + " bonds\n")
    data_file.write(str(n_angles) + " angles\n")
    data_file.write("0 dihedrals\n")
    data_file.write("0 impropers\n")
    data_file.write("\n")

    if version == 5:
        if HPS:
            data_file.write("23 atom types\n")
        else:
            data_file.write("43 atom types\n")

        if ENM:
            data_file.write(str(t - 1) + " bond types\n")
        else:
            data_file.write("3 bond types\n")

    else:

        if HPS:
            data_file.write("22 atom types\n")
        else:
            data_file.write("42 atom types\n")

        if ENM:
            data_file.write(str(t - 1) + " bond types\n")
        else:
            data_file.write("3 bond types\n")
    data_file.write("0 angle types\n")
    data_file.write("\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " xlo xhi\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " ylo yhi\n")
    data_file.write("" + str(-max_r) + " " + str(max_r) + " zlo zhi\n")

    bp_file = open("DNA_sequence.txt", "w")
    bp_file.write("# " + str(n_dna) + "\n")
    data_file.write("\n")
    data_file.write("Atoms\n")
    data_file.write("\n")
    i = 1
    m = 0

    if version == 1 or version == 3 or version == 4 or version == 5:

        mol_starts = list()
        mol_ends = list()
        for alt_proteins in nucls:
            p = 1
            for protein in alt_proteins:
                if p == 1:
                    m = m + 1
                else:
                    m = m

                mol = m
                s = 1
                for subsection in protein:
                    start_i = i

                    for atom in subsection:
                        a_type = index_dict_rev[atom[0]]

                        if s == 2:
                            # rigid region
                            if not HPS:
                                a_type = a_type + 20

                        charge = charge_dict[atom[0]]
                        # atom id, atom type, x, y, z ellipsoid flag, density,
                        # molecule, charge
                        data_file.write(
                            str(i) + " " + str(a_type) + " " + str(atom[2][0]) + " " +
                            str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                            str(mol) + " " + str(charge) + "\n")

                        i = i + 1
                    end_i = i - 1

                    if s == 2:
                        print(str(start_i) + ":" + str(end_i))
                        mol_starts.append(start_i)
                        mol_ends.append(end_i)
                    s = s + 1
                if LH:
                    if p == 1:
                        m = m + 1
                p = p + 1

            print("\n")
        m = m
        mol = m
        if version == 5 and LH:
            mmm = 1
            for atom in dummy_atom:
                if HPS:
                    X_TYPE = 23
                else:
                    X_TYPE = 43
                data_file.write(
                    str(i) + " " + str(X_TYPE) + " " + str(atom[2][0]) + " " +
                    str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                    str(mmm) + " " + str(0.0) + "\n")
                i = i + 1
                mmm = mmm + 1

        ellipsoid_list = list()
        first_dna_mol = mol + 1
        for bead in dna_full:
            # atom id, atom type, x, y, z ellipsoid flag, density, molecule, charge
            if bead[0] == "D":
                a_type = D_type
                el = 1
                charge = 0.0
                ellipsoid_list.append(i)
                bp_file.write(str(i) + " " + bead[1][0] + bead[1][2] + "\n")
                mol = mol + 1
            else:
                a_type = P_type
                el = 0
                charge = -1.0

            data_file.write(
                str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                str(mol) + " " + str(charge) + "\n")

            i = i + 1
        last_dna_mol = mol
    else:  # need to group entire nucleosome protein + dna into one molecule
        m = 1
        mol_starts = list()
        mol_ends = list()

        for alt_proteins in nucls:
            p = 1
            for protein in alt_proteins:
                # if p == 1:
                #     m = m + 1
                # else:
                #     m = m

                mol = m
                s = 1
                for subsection in protein:
                    start_i = i

                    for atom in subsection:
                        a_type = index_dict_rev[atom[0]]

                        if s == 2:
                            # rigid region
                            if not HPS:
                                a_type = a_type + 20

                        charge = charge_dict[atom[0]]
                        # atom id, atom type, x, y, z ellipsoid flag, density,
                        # molecule, charge
                        data_file.write(
                            str(i) + " " + str(a_type) + " " + str(atom[2][0]) + " " +
                            str(atom[2][1]) + " " + str(atom[2][2]) + " 0  1 " +
                            str(mol) + " " + str(charge) + "\n")

                        i = i + 1
                    end_i = i - 1

                    if s == 2:
                        print(str(start_i) + ":" + str(end_i))
                        mol_starts.append(start_i)
                        mol_ends.append(end_i)
                    s = s + 1
                # if p == 1:
                #     m = m + 1
                p = p + 1

            m = m + 1

            print("\n")
        m = m
        mol = m
        ellipsoid_list = list()

        dna_nucl_sections = list()
        dna_nucl_atom_ids = list()
        nrl = int(dna_len_inital / 3)

        mid_point = nrl / 2

        start_point = int(mid_point - 147 / 2)
        end_point = int(mid_point + 147 / 2)

        for N in range(NUM_NUCLS):
            dna_nucl_sections.append([start_point * 3 + N * nrl * 3, end_point * 3 + N * nrl * 3 + 2])

            dna_nucl_atom_ids.append([dna_nucl_sections[-1][0] + i, dna_nucl_sections[-1][1] + i])
        j = 0
        for bead in dna_full:
            # atom id, atom type, x, y, z ellipsoid flag, density, molecule, charge
            if bead[0] == "D":
                a_type = D_type
                el = 1
                charge = 0.0
                ellipsoid_list.append(i)
                bp_file.write(str(i) + " " + bead[1][0] + bead[1][2] + "\n")
                mol = mol + 1
            else:
                a_type = P_type
                el = 0
                charge = -1.0
            n = 1

            written = False
            for nucl_dna in dna_nucl_sections:
                if (j <= nucl_dna[1] and j >= nucl_dna[0]):
                    m = n

                    data_file.write(
                        str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                        str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                        str(m) + " " + str(charge) + "\n")

                    written = True

                n = n + 1

            if not written:
                data_file.write(
                    str(i) + " " + str(a_type) + " " + str(bead[2][0]) + " " +
                    str(bead[2][1]) + " " + str(bead[2][2]) + " " + str(el) + " 1 " +
                    str(mol) + " " + str(charge) + "\n")

            i = i + 1
            j = j + 1

    bp_file.close()

    # if seq_type == 5:
    #     seq_file = open("DNA_sequence.txt",'r')
    #     data = list()
    #     for line in seq_file:
    #         data.append(line)
    #
    #     seq_file.close()
    #     s_len = len(data)
    #
    #     data_array = list()
    #     for k in range(1,s_len):
    #         split_s = data[k].split()
    #
    #         data_array.append([split_s[0], split_s[1]])
    #
    #     data_array = np.array(data_array)
    #     diff = new_nrl_length - 211
    #     ncp147 = data_array[int(np.floor(diff/2)):-int(np.ceil(diff/2)),1]
    #
    #
    #
    #     v1 = data_array[:]
    #     v1[:211,1] = ncp147[:]
    #     v1[-211:,1] = ncp147[:]
    #     v1[211:-211] = "AT"
    #
    #     v2 = v1[:]
    #
    #     v2[211:-211] = "TA"
    #
    #     v1_file = open("DNA_sequence_2.txt",'w')
    #     v2_file = open("DNA_sequence_3.txt", 'w')
    #     v1_file.write(data[0])
    #     v2_file.write(data[0])
    #
    #     for k in range(0,s_len-1):
    #         v1_file.write(v1[k,0] + " " + v1[k,1]+"\n")
    #         v2_file.write(v2[k,0] + " " + v2[k,1]+"\n")
    #     v1_file.close()
    #     v2_file.close()

    if seq_type == 5:
        seq_file = open("DNA_sequence.txt", 'r')
        data = list()
        for line in seq_file:
            data.append(line)

        seq_file.close()
        s_len = len(data)

        col_1 = list()
        col_2 = list()

        for k in range(1, s_len):
            split_s = data[k].split()

            col_1.append(split_s[0])
            col_2.append(split_s[1])

        v1_file = open("DNA_sequence_2.txt", 'w')
        v1_file.write(data[0])

        new_col_2 = list()
        for line in col_2:
            new_col_2.append(line)

        for i in range(left_l, left_l + 211):
            new_col_2[i] = col_2[-(211 + left_l) + (i - left_l)]
            new_col_2[-(211 + left_l) + (i - left_l)] = col_2[i]

        for k in range(0, s_len - 1):
            v1_file.write(col_1[k] + " " + new_col_2[k] + "\n")
        v1_file.close()
    
    if not LH:
        nucl_string = "group nucl id"
        for start, end in zip(mol_starts, mol_ends):
            nucl_string = nucl_string + " " + str(start) + ":" + str(end)
    else:
        nucl_string = "group nucl id"
        for n in range(NUM_NUCLS):
            for start, end in zip(mol_starts[(n*9)+1:((n+1)*9)], mol_ends[(n*9)+1:((n+1)*9)]):
                nucl_string = nucl_string + " " + str(start) + ":" + str(end)
        nucl_string=nucl_string+"\ngroup lh id"
        for n in range(NUM_NUCLS):
            start = mol_starts[(n*9)]
            end   = mol_ends[(n*9)]
            nucl_string = nucl_string + " " + str(start) + ":" + str(end)
 
    if version == 2:
        print("here")
        for dna_bloc in dna_nucl_atom_ids:
            nucl_string = nucl_string + " " + str(dna_bloc[0]) + ":" + str(dna_bloc[1])

    print(nucl_string, file=in_file)

    data_file.write("\n")
    data_file.write("Ellipsoids\n")
    data_file.write("\n")
    # print("group dna_all id " + str(ellipsoid_list[0]) + ":" + str(ellipsoid_list[-1]), file=in_file)
    # print("group notdna subtract all dna_all", file=in_file)
    # print("group tails  subtract notdna nucl", file=in_file)
    # print("group dna id " + str(ellipsoid_list[0]) + ":" + str(ellipsoid_list[-1]) + ":3", file=in_file)
    a = 0
    for bead in dna_full:
        if bead[0] == "D":
            data_file.write(
                str(ellipsoid_list[a]) + " " + str(RX) + " " + str(RY) + " " +
                str(RZ) + " " + str(bead[3][0]) + " " + str(bead[3][1]) + " " +
                str(bead[3][2]) + " " + str(bead[3][3]) + "\n")
            a = a + 1

    data_file.write("\n")
    data_file.write("Bonds\n")
    data_file.write("\n")
    for bond in bonds:
        data_file.write(
            str(bond[0]) + " " + str(bond[1]) + " " + str(bond[2]) + " " +
            str(bond[3]) + "\n")

    data_file.write("\n")
    # data_file.write("Angles\n")
    # data_file.write("\n")
    # for angle in angles:
    #     data_file.write(
    #         str(angle[0]) + " " + str(angle[1]) + " " + str(angle[2]) + " " +
    #         str(angle[3]) + " " + str(angle[4]) + "\n")
    # data_file.close()

    # print("group end1 id " + str(ellipsoid_list[0]), file=in_file)
    # print("group end2 id " + str(ellipsoid_list[-1]), file=in_file)
    # print("group dna_P subtract dna_all dna", file=in_file)

    if HPS:
        print("group dna type 21", file=in_file)
        print("group ps type 22", file=in_file)
    else:
        print("group dna type 41", file=in_file)
        print("group ps type 42", file=in_file)

    if version == 1 or version == 3 or version == 4 or version == 5:
        print("group dna_all union dna ps", file=in_file)
        print("group notdna subtract all dna_all", file=in_file)
        print("group tails subtract notdna nucl", file=in_file)
        print("neigh_modify exclude molecule/intra nucl", file=in_file)
        if LH:
            print("neigh_modify exclude molecule/intra lh", file=in_file)
        print("neigh_modify exclude molecule/intra dna_all ", file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(P_type), file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(D_type), file=in_file)

        first_dna_id = ellipsoid_list[0]
        last_dna_id = ellipsoid_list[-1]

        # nrl = (last_dna_id - first_dna_id)/3
        # diff = nrl - 211
        # half = int(diff/2)*3
        half = 0
        print('variable ex equal "sqrt((x[' + str(first_dna_id + half) + '] - x[' + str(
            last_dna_id - half) + '])^2 + (y[' + str(first_dna_id + half) + '] - y[' + str(
            last_dna_id - half) + '])^2 + (z[' + str(first_dna_id + half) + '] - z[' + str(
            last_dna_id - half) + '])^2)"', file=in_file)
        # print('group end1 molecule '+str(first_dna_mol),file=in_file)
        # print('group end2 molecule '+str(last_dna_mol),file=in_file)
        print('group end1 id ' + str(first_dna_id), file=in_file)
        print('group end2 id ' + str(last_dna_id), file=in_file)
        print('group P1 id ' + str(first_dna_id + 1), file=in_file)
        print('group P2 id ' + str(first_dna_id + 2), file=in_file)
        print('group p3 id ' + str(last_dna_id + 1), file=in_file)
        print('group p4 id ' + str(last_dna_id + 2), file=in_file)
        # compute the dna gyration
        dnagyr = get_gyration(dna_full)
        print(dnagyr)
        
        
        # make colvars file
        #print("mkaing colvars file")
        #colvarfile=open("colvars_test.txt","w")
        #colvarfile.write("colvar {\n")
        #colvarfile.write("name gyr\n")
        #colvarfile.write("gyration {\n")
        #colvarfile.write("atoms {\n")
        #colvarfile.write("atomNumbersRange 1-" + str(last_dna_id)+"\n")
        #colvarfile.write("}\n")
        #colvarfile.write("}\n")
        #colvarfile.write("lowerBoundary 1\n")
        #colvarfile.write("upperBoundary "+str(dnagyr)+"\n")
        ##colvarfile.write("upperBoundary 1000\n")
        #colvarfile.write("}\n")
        
        
        
        #L=len(mol_starts)
        #nucl_id_lists=list()
        #for n in range(NUM_NUCLS):
        #    if LH:
        #        tem1 = mol_starts[1+n*9:(n+1)*9]
        #        tem2 = mol_ends[1+n*9:(n+1)*9]
        #    else:
        #        tem1 = mol_starts[n*8:(n+1)*8]
        #        tem2 = mol_ends[n*8:(n+1)*8]
        #    nucl_id_lists.append([tem1,tem2])
        #print(nucl_id_lists)
 
        #
        ## all distance pairs
        ## for N nucleosomes we have 
        #a=1
        #for n in range(NUM_NUCLS):
        #    for m in range(n+1,NUM_NUCLS):     
        #        colvarfile.write("colvar {\n")
        #        colvarfile.write("name d"+str(a)+"\n")
        #        a=a+1
        #        colvarfile.write("distance {\n")
        #        b=1
        #        for g in [n,m]:
        #            # the 2 groups
        #            colvarfile.write("group"+str(b)+" {\n")
        #            b=b+1
        #            starts=nucl_id_lists[g][0]
        #            ends=nucl_id_lists[g][1]
        #            for start, end in zip(starts,ends):
        #                rstring ="atomNumbersRange " + str(start) + "-" + str(end)+"\n"
        #                colvarfile.write(rstring)
        #            colvarfile.write("}\n")
        #        colvarfile.write("}\n")
        #        colvarfile.write("lowerBoundary 50\n")
        #        colvarfile.write("upperBoundary 500\n")
        #        colvarfile.write("}\n")

     
        
        
        
        
        
        
        
        ## do the angles, for N nucls we have N-2 angles
        #L=len(mol_starts)
        #nucl_id_lists=list()
        #for n in range(NUM_NUCLS):
        #    if LH:
        #        tem1 = mol_starts[1+n*9:(n+1)*9]
        #        tem2 = mol_ends[1+n*9:(n+1)*9]
        #    else:
        #        tem1 = mol_starts[n*8:(n+1)*8]
        #        tem2 = mol_ends[n*8:(n+1)*8]
        #    nucl_id_lists.append([tem1,tem2])
        #print(nucl_id_lists)
        #for n in range(0,NUM_NUCLS-2):
        #    colvarfile.write("colvar {\n")
        #    colvarfile.write("name a"+str(n+1)+"\n")
        #    colvarfile.write("angle {\n")
        #    for g in range(3):
        #        # the three groups
        #        colvarfile.write("group"+str(g+1)+" {\n")
        #        starts=nucl_id_lists[g+n][0]
        #        ends=nucl_id_lists[g+n][1]
        #        for start, end in zip(starts,ends):
        #            rstring ="atomNumbersRange " + str(start) + "-" + str(end)+"\n"
        #            colvarfile.write(rstring)
        #        colvarfile.write("}\n")
        #    colvarfile.write("}\n")
        #    colvarfile.write("lowerBoundary 0\n")
        #    colvarfile.write("upperBoundary 180\n")
        #    colvarfile.write("}\n")

        ##metadynamics {
        ##name meta
        ##colvars a1 a2 a3 a4
        ##hillWeight 100
        ##hillWidth 1
        ##}


        #colvarfile.write("harmonic {\n")
        #colvarfile.write("name harm\n")
        #colvarfile.write("colvars gyr\n")
        #colvarfile.write("centers 10\n")
        #colvarfile.write("forceConstant 0.05\n")
        #colvarfile.write("}\n")
        #
        #
        #
        #colvarfile.close()
        
        
        
        
        #if not RESTART:
        #    settings = open("lmp_run_settings_nucl.txt", 'r')
        #    for line in settings:
        #        if not version == 5:
        #            if line == "comm_modify cutoff 80\n":
        #                line = "\n"

        #        in_file.write(line)
        #else:
        #    settings = open("lmp_run_settings_nucl_restart.txt", 'r')
        #    for line in settings:
        #        if not version == 5:
        #            if line == "comm_modify cutoff 80\n":
        #                line = "\n"
        #        in_file.write(line)
    else:
        print("group dna_all_temp union dna ps", file=in_file)
        print("group dna_all subtract dna_all_temp nucl", file=in_file)
        print("group temp subtract all nucl", file=in_file)
        print("group tails subtract temp dna_all", file=in_file)
        print("group dna_only intersect dna_all dna", file=in_file)

        print("neigh_modify exclude molecule/intra nucl", file=in_file)
        print("neigh_modify exclude molecule/intra dna_all ", file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(P_type), file=in_file)
        print("neigh_modify exclude type " + str(D_type) + " " + str(D_type), file=in_file)

        #settings = open("lmp_run_settings_v2.txt", 'r')
        #for line in settings:
        #    in_file.write(line)

    in_file.close()

def get_gyration(dna_full):
    pos=list()
    for bead in dna_full:
        if bead[0] == "D":
            pos.append(bead[2])
    num_d=len(pos)
    pos=np.array(pos)
    rcom=np.mean(pos,axis=0)
    #print("RCOM = ",rcom)
    R=0.0
    for p in pos:
        d=np.sum((p-rcom)**2)
        #print(d)
        R=R+d
    return np.sqrt(R/(num_d))



def update_from_dump(alt_proteins, dna_full, n_dna, dump_name,ids):
    lines = [line.rstrip('\n') for line in open(dump_name)]
    print("read in dump")
    num_atoms=int(lines[3])
    #print(dna_full)
    coords =list()
    for i in range(9,9+num_atoms):
        coords.append(np.array(lines[i].split(),dtype=float))
    #print(ids)
    protein_list = list()
    Ds = list()
    Ps = list()
    for coord in coords:
        btype = int(coord[-1])
        #print(btype)
        if btype == ids[0]:
            Ds.append(coord[1:8])
        elif btype == ids[1]:
            Ps.append(coord[1:4])
        else:
            protein_list.append(coord[1:4])
    
    #print(Ds)
    for i in range(0,len(Ds)):
        #print(dna_full[i*3])
        #print(Ds[i])

        dna_full[i*3][2]    = Ds[i][:3]
        dna_full[i*3][3]    = Ds[i][3:]
        
        #print(dna_full[i*3])
        #print("")
        dna_full[i*3+1][2]  = Ps[i*2][:]
        dna_full[i*3+2][2]  = Ps[i*2+1][:]
    i=0
    for protein in alt_proteins:
        for subsection in protein:
            for bead in subsection:
                bead[2] = protein_list[i][:]
                i=i+1

    return alt_proteins, dna_full, n_dna

def create_bead_string(LH,num_nucls,nrl,rigid,seq,ishremd, path_to_lammps):
    """ Make a single nucleosome
    
    Arguments:
        LH: bool - inlcude linker histone (LH) <True> or not <False>
        num_nucls: int - number of nucleosomes
        nrl: int - Nucleosome repeat length, needs to >147, may give unexpected
                    results when > 211.
        rigid: string - <"full"> for fully bonded DNA to nucl,
                        <"none"> for non DNA-nucl bonding
         seq: string  - <"ncp147"> for squence of ncp147
                        <"AAAA"> for poly A
         ishremd: bool - do the setup for HREMD <True> or not <False>
         path_to_lammps: string - path to lammps exe
    """

    # 1. single nucleosome, verions 4, pull part to make vertical
    create_nucl(0.1, 0.01, False, 1, 4, True, 1, False, True, 0,0)
    time.sleep(5)
    subprocess.call('cat in.nucl lmp_run_settings_nucl_bead_string_setup.txt > in.setup',shell=True)
    time.sleep(5)
    subprocess.call('mpirun -np 8 '+ path_to_lammps+' -in in.setup',shell=True)
    print("run setup \n")
    time.sleep(5)
    subprocess.call('tail -n $(($(head -n 10 bead_string_setup.dump | grep -A1 "NUMBER OF ATOMS" | tail -1) + 9))  bead_string_setup.dump > bead_string_setup_last_frame.dump', shell=True)
    print("made dump\n")
    time.sleep(20)
    
    # now we use it as the basis for the generation of the full chromatin
    lower=-int(np.ceil((211-nrl)/2))
    upper=-int(np.floor((211-nrl)/2))
    if rigid=="full" or rigid=="partial":
        version=3
        if rigid=="full":
            RESTN=147
        else:
            RESTN=74
    else:
        version=1
        RESTN=0
    if seq == "ncp147":
        seqint = 1


    if seq == "AAAA":
        seqint = 4
    # create the rigid version
    #create_nucl_from_dump(0.1, 0.01, False, num_nucls, 2, LH, seqint, False, True, upper,lower)
    
    # create the normal version
    create_nucl_from_dump(0.1, 0.01, False, num_nucls, version, LH, seqint, False, True, upper,lower,RESTN,ishremd)
    #create_nucl_from_dump(0.1, 0.01, False, num_nucls, version, LH, seqint, True, True, upper,lower)

    time.sleep(5)
    if ishremd:
        subprocess.call('cat hremd_replica_settings.txt in.nucl hremd_run_settings.txt > in.hremd', shell=True)
        subprocess.call("sed -i 's/variable a equal 8.0//g' in.hremd",shell=True)
        
    else:
        subprocess.call('cat in.nucl run_settings.txt > in.run', shell=True)
        subprocess.call('mv in.run in.nucl',shell=True)




# def create_normal(LH,num_nucls,nrl,rigid,seq):
#     lower = -int(np.ceil((211 - nrl) / 2))
#     upper = -int(np.floor((211 - nrl) / 2))

#     if rigid:
#         version = 3

#     else:
#         #if LH:
#         #    version = 5
#         #else:
#         version = 1

#     if seq == "ncp147":
#         seqint=1

#     if seq == "AAAA":
#         seqint=4

#     create_nucl(0.1, 0.01, False, num_nucls, version, LH, seqint, False, True, upper, lower)
#     create_nucl(0.1, 0.01, False, num_nucls, version, LH, seqint, True, True, upper, lower)


def create_single(LH,nrl,rigid,seq):
    """ Make a single nucleosome
    
    Arguments:
        LH: bool - inlcude linker histone (LH) <True> or not <False>
        nrl: int - Nucleosome repeat length, needs to >147, may give unexpected
                    results when > 211.
        rigid: string - <"full"> for fully bonded DNA to nucl,
                        <"none"> for non DNA-nucl bonding
         seq: string  - <"ncp147"> for squence of ncp147
                        <"AAAA"> for poly A
    """
    
    num_nucls=1
    lower = -int(np.ceil((211 - nrl) / 2))
    upper = -int(np.floor((211 - nrl) / 2))

    if rigid == "full":
        version = 3

    else:
        #if LH:
        #    version = 5
        #else:
        version = 1

    if seq == "ncp147":
        seqint=1

    if seq == "AAAA":
        seqint=4

    create_nucl(0.1, 0.01, False, num_nucls, version, LH, seqint, False, True, upper, lower)
    # create_nucl(0.1, 0.01, False, num_nucls, version, LH, seqint, True, True, upper, lower)






