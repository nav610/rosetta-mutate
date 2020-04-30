from pyrosetta import *
from pyrosetta.toolbox import mutate_residue
init()

aa_long = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
aa_short = {'ALA' : 'A', 'ARG' : 'R', 'ASN' : 'N', 'ASP' : 'D', 'CYS' : 'C', 'GLN' : 'Q', 'GLU' : 'E', 'GLY' : 'G', 'HIS' : 'H','HIS_D' : 'H', 'ILE' : 'I', 'LEU' : 'L', 'LYS' : 'K', 'MET' : 'M', 'PHE' : 'F', 'PRO' : 'P', 'SER' : 'S', 'THR' : 'T', 'TRP' : 'W', 'TYR' : 'Y', 'VAL' : 'V'}
chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

f = open("surface.txt",'r')
res_list = [int(i) for i in f.readlines()]
f.close()

f = open("energies.txt",'w')

pose_org = pose_from_pdb("nlrp6_symm.pdb")
scorefx = get_fa_scorefxn()
f.write("org {}\n".format(scorefx(pose_org)))
errorlist = []

#loop starts here; first layer removed, would include looping through res_list
for res in res_list: 
    try:
        cur_long = pose_org.residue(res-13).name()
        cur_short = aa_short[cur_long]
        #for mut in aa_long: #["ALA"]
        for mut in aa_long:
            if mut != cur_long:
                new_pose = Pose()
                new_pose.assign(pose_org)
                for c in chains:
                    num= new_pose.pdb_info().pdb2pose(c,res)
                    mutate_residue(new_pose,num, aa_short[mut],5,scorefx) 
                f.write("{}{}{}.pdb {}\n".format(cur_short,res,aa_short[mut],scorefx(new_pose)))
                new_pose.dump_pdb("full_results/{}{}{}.pdb".format(cur_short,res,aa_short[mut]))
    except:
        errorlist.append(res)
        continue
f.close()
for i in errorlist:
    print("Res {} did not finish".format(i))
