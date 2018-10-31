import ROOT as rt
from time import sleep
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle

from ClusterMatcher import *

#Main function:

def TrackMatch(file, dR, MomentumThreshold, HadronsNotQuarks=False, Plot=False, Axes=None, Save=False, dR_dist = False, EarlyBreak=0, manual_selection = False):
	"""returns unique ID, vertices and 3-momenta of all tracks that lie inside the dR-cone of a b-particle trajectory; optionally it returns also a 3D-plot
		

	Inputs:
		file: 			full root TFile
		dR: 			Delta R region around particle trajectory in which clusters should count as hit
		Momentumthreshold:	momentum threshold for b-particles to be counted
		HadronsNotQuarks:       if set as True, the algorithm will focus on B-hadrons instead of b-quarks
		Plot:			if set as True, the function will return a 3D-plot
		Axes:                   pre-initialized 3D-plot axes. Only necessary if Plot==True
		Save:			if set as True, the function will save the data to a .pkl file for later use
		dR_dist:		if set as True, the function will return a histrogram showing the distribution of delta R between pixel clusters and the corresponding trajectory
		EarlyBreak:             non zero integer which denotes the number of events after which the algorithm should stop
		manual_selection:	optional input: tuple (nEvent,nParticle) for selecting a specific event
	
	Outputs:
		list of tuples where each contains a tuple with a unique identification followed by the cartesian coordinates of the vertex and then the momentum vector: ((nEvent,nParticle,nTrack),(vx_x,vx_y,vx_z),(px,py,pz))
		"""
	#colorstring for plots:
	#color = ['b','g','r','c','m','y','k','w']
	hsv = plt.get_cmap('hsv')
	color = hsv(np.linspace(0,1.0,12))
	c = 0 #initializing color index
	res = 50 #number of points in trajectory

	# open tree file
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	MatchedTracks = []
	#hist = rt.TH1D('B-hadron Pt','B-hadron Pt',40,0.,510.)
		
	for i in xrange(N):
    		if i % 50 == 0: print "Working on event " ,i
		if EarlyBreak >0 and i >= EarlyBreak: break
		if manual_selection != False:
			if i < manual_selection[0]: 
				continue
			elif i > manual_selection[0]:
				break
    		tree.GetEntry(i)

		for j in range(0,tree.nJets):
        		jVector = rt.TLorentzVector()
        		jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
        		for k in range(0,tree.nGenParticles):
				if manual_selection != False:
					if k < manual_selection[1]: 
						continue
					elif k > manual_selection[1]:
						break

				if HadronsNotQuarks == False:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) == 5
					statusCriterion = tree.genParticle_status[k] == 23
				else:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600
					statusCriterion = tree.genParticle_status[k] == 2 
            			if statusCriterion and pdgCriterion:
                    			pVector = rt.TLorentzVector()
                    			pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                        			tree.genParticle_phi[k],tree.genParticle_mass[k])
                    			delR = jVector.DeltaR(pVector)
					#if delR < 0.3:
					#	hist.Fill(tree.genParticle_pt[k])                       
                    			if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
                        			v_p = normalize(np.array([pVector[0]/tree.genParticle_mass[k], pVector[1]/tree.genParticle_mass[k], \
                                   			pVector[2]/tree.genParticle_mass[k]]))
                        			phi = PolarPhi(v_p[0],v_p[1])
						theta = Theta(v_p[0],v_p[1],v_p[2])
						eta = Eta(theta)
						if np.sqrt(tree.genParticle_decayvx_x[k]**2 + tree.genParticle_decayvx_y[k]**2) > 4:
							print "happens at event nr {}, particle nr {}".format(i, k)

						if Plot == True:    
							t_max = TrajectoryLength(theta,v_p)
							PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,Axes,t_max,res,color[c],1,'--')
                                               	
                        			NearTracks = [] 
						#x3 = 0
						for nTrack in xrange(0,tree.nTracks): #finding all tracks inside deltaR<dR
							VertexTheta,VertexPhi = ShiftedThetaPhi(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack],tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
							DRT = DeltaR_eta(eta,tree.track_eta[nTrack],phi,tree.track_phi[nTrack])
							dxyz = np.sqrt((tree.track_vx_x[nTrack]-tree.genParticle_vx_x[k])**2 + (tree.track_vx_y[nTrack]-tree.genParticle_vx_y[k])**2 + (tree.track_vx_z[nTrack]-tree.genParticle_vx_z[k])**2)
							if DRT<dR: # and dxyz<0.5:
								#r = np.sqrt(tree.track_vx_x[nTrack]**2 + tree.track_vx_y[nTrack]**2)
								#if r > 4.5:
								#	x2 += 1
								#	x3 += 1
								#track_vx_r_hpt.Fill(r)
								#dxyz_hist.Fill(dxyz)
								NearTracks.append(((i,k,nTrack),(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack]),(tree.track_px[nTrack],tree.track_py[nTrack],tree.track_pz[nTrack])))
						if  NearTracks == []:
							break	
						else:
							MatchedTracks.append(NearTracks) 
                        				if Plot == True:
								
								for entry in NearTracks:
									PlotTrajectory(entry[1],entry[2],Axes,0.7*t_max,res,color[c],0.5,'-')
								Axes.scatter(tree.genParticle_decayvx_x[k],tree.genParticle_decayvx_y[k],tree.genParticle_decayvx_z[k], color='black',s=25)
						if Plot == True:
							if c != len(color)-1:
                        					c += 1
                        				else:
                            					c = 0
	
	if HadronsNotQuarks == False:
		particleString = 'b-quarks'
	else:
		particleString = 'B-hadrons'
	print "Total Number of high pt "+particleString+": ", len(MatchedTracks)
	nMatchedTracks = 0
	for entry in MatchedTracks:
		nMatchedTracks += len(entry)
	print "Total number of matched tracks ", nMatchedTracks
	if Plot == True:
		plt.savefig("Trajectories_Tracks.png")
		plt.show()
	if Save == True:
		with open("MatchedTracksDR"+str(dR)+"on"+str(particleString)+".pkl", 'w') as f:
			pickle.dump(MatchedTracks, f)

	return MatchedTracks

if __name__ == '__main__':
	
	file = rt.TFile("flatTuple.root",'READ')

	dR = 0.1 #DeltaR threshold for counting tracks
	MomentumThreshold = 350

	with open("Grid.pkl",) as f:    #open coordinates of DetUnits for visual reference
                Grid = pickle.load(f)

	ax = Initialize3DPlot('Particle Trajectories', 'x', 'y', 'z', grid=Grid)
	MatchedTracks =  TrackMatch(file, dR, MomentumThreshold, HadronsNotQuarks=True, Plot=True, Axes=ax, Save=False, dR_dist = False, EarlyBreak=10)

	#with open("MatchedTracksDR0.1onB-hadrons.pkl",) as f:	#take data from file instead of running entire function
	#	MatchedTracks = pickle.load(f)


	
	'''
	Signal_4TeV_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_noPU/180628_092048/0000/flatTuple_{}.root'
	file = rt.TFile.Open(Signal_4TeV_noPU_String.format(20))
	tree = file.Get("demo/tree")
	N = tree.GetEntries()
	ax = Initialize3DPlot('B-Hadron Evolution', 'x', 'y', 'z', grid=None )	
	tree.GetEntry(0)
	#particle = tree.genParticle[0]
	#import pdb; pdb.set_trace()


	def PlotTracks(nEvent,nParticle,color):
		i = nEvent
		k = nParticle
		tree.GetEntry(i)
		#print "Pt =", tree.genParticle_pt[k]
		pVector = rt.TLorentzVector()
        	pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k],tree.genParticle_phi[k],tree.genParticle_mass[k])
		v_p = normalize(np.array([pVector[0]/tree.genParticle_mass[k], pVector[1]/tree.genParticle_mass[k], \
        		pVector[2]/tree.genParticle_mass[k]]))
		phi = PolarPhi(v_p[0],v_p[1])
		theta = Theta(v_p[0],v_p[1],v_p[2])
		eta = Eta(theta)
		t_max = TrajectoryLength(theta,v_p)
		PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,ax,0.8*t_max,50,'green',1,'--')
		ax.scatter(tree.genParticle_decayvx_x[k],tree.genParticle_decayvx_y[k],tree.genParticle_decayvx_z[k],c='darkolivegreen',s=20,linewidths=0.1)
		print "Event", nEvent, ",Particle", nParticle
		print "B-hadron vertex at", tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]
		print "B-hadron decay-vertex at", tree.genParticle_decayvx_x[k],tree.genParticle_decayvx_y[k],tree.genParticle_decayvx_z[k]
		NearTracks = [] 
		for nTrack in xrange(0,tree.nTracks): #finding all tracks inside deltaR<dR
			VertexTheta,VertexPhi = ShiftedThetaPhi(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack],tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
			DRT = DeltaR_eta(eta,tree.track_eta[nTrack],phi,tree.track_phi[nTrack])
			if DRT<dR:
				NearTracks.append(((i,k,nTrack),(tree.track_vx_x[nTrack],tree.track_vx_y[nTrack],tree.track_vx_z[nTrack]),(tree.track_px[nTrack],tree.track_py[nTrack],tree.track_pz[nTrack])))
		for nT, entry in enumerate(NearTracks):
			print "Track", nT, ", vertex at", entry[1][0], entry[1][1], entry[1][2]
			PlotTrajectory(entry[1],entry[2],ax,0.7*t_max,50,color,0.5,'-')
		#plt.show()
	
	ClustersX = [-2.1151912212371826, -1.9214503765106201, -3.309929609298706, -3.2997121810913086, -3.297450542449951, -3.3190650939941406, -3.3190648555755615, -3.008465051651001, -5.136828422546387, -4.964303970336914, -4.946408271789551, -4.9160847663879395, -4.912755966186523, -4.901236534118652, -4.467060565948486]
	ClustersY = [4.139054298400879, 4.209573268890381, 6.241778373718262, 6.246328353881836, 6.247335910797119, 6.237710952758789, 6.237710475921631, 6.3760175704956055, 9.06138801574707, 9.144474029541016, 9.153093338012695,9.16769790649414, 9.169301986694336, 9.174849510192871, 9.383955955505371]
	ClustersZ = [9.683456420898438, 9.266300201416016, 11.82179069519043, 11.784566879272461, 11.762925148010254, 11.799983978271484, 11.82470703125, 11.171915054321289, 14.886789321899414, 14.806136131286621, 14.735239028930664, 14.6824369430542, 14.646608352661133, 14.65557861328125, 13.795677185058594]                                                                                                                                         
	PlotTracks(0,5,'red')
	print "particle pddId = ", tree.genParticle_pdgId[5]
	ax.scatter(ClustersX[:2], ClustersY[:2], ClustersZ[:2], c='blue',s=7,linewidths=0.1)
	ax.scatter(ClustersX[2:8], ClustersY[2:8], ClustersZ[2:8], c='blue',s=7,linewidths=0.1)
	ax.scatter(ClustersX[8:], ClustersY[8:], ClustersZ[8:], c='blue',s=7,linewidths=0.1)

	ax.set_xlim3d(-10,1)
	ax.set_ylim3d(-1.15)
	ax.set_zlim3d(0,25)
	plt.show()
	'''
	

	'''27, 6: 27, 9: 195, 5'''

	'''Pt = 620.925048828, jet_bTag = 0.845161676407, Hits: 2,6,7'''
	














