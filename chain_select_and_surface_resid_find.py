from pyrosetta import *
import argparse
init()

parser = argparse.ArgumentParser() 
parser.add_argument("-pdb","--pdb",help = "pdb") 
args = parser.parse_args() 

sasa=pyrosetta.rosetta.core.scoring.sasa.SasaCalc()

pose = pose_from_pdb(args.pdb)
sasa.calculate(pose)
sasa_residue = list(sasa.get_residue_sasa())

f = open("newBfactor.txt",'w')
for i in range(len(sasa_residue)):
	f.write(str(sasa_residue[i]))
	f.write("\n")

f = open("newBfactor_step.txt",'w')
s = open("surface.txt",'w')
threshold=float(50.0)
for i in range(len(sasa_residue)):
	if float(sasa_residue[i])>threshold:
		f.write(str(1) + "\n")
		s.write(str(pose.pdb_info().pose2pdb(i+1).split( )[0])+"\n") #save as pdb format pose2pdb



	else: f.write(str(0)+ "\n")


