from DisHist import *
import numpy as np
import ROOT as rt
from time import sleep


def Histogramize(Histograms, ran, title, xlabel, ylabel, Save=False,Normalize=True, t_sleep=0):
	canvas = rt.TCanvas('canvas','canvas',600,600)
	canvas.SetTitle(title)
	if len(Histograms) > 1: 
		rt.gStyle.SetOptStat(0)
		legend = rt.TLegend(0.9,0.9,0.65,0.75)
	for nr, Histogram in enumerate(Histograms):
		Histogram[0].SetLineColor(nr+2)
		if nr == 0:
			Histogram[0].GetXaxis().SetTitle(xlabel)
			Histogram[0].GetYaxis().SetTitle(ylabel)
			Histogram[0].GetYaxis().SetTitleOffset(1.5)
			if Normalize:
				Histogram[0].DrawNormalized()
			else:
				Histogram[0].Draw()
		else:
			if Normalize:
				Histogram[0].DrawNormalized("SAME")
			else:
				Histogram[0].Draw("SAME")
		if len(Histograms)>1:
			legend.AddEntry(Histogram[0],Histogram[1])
	if len(Histograms)>1: legend.Draw()
	if Save: canvas.SaveAs(title+".png")
	sleep(t_sleep)

def LoadHistogram(Histogram, filename):
	with open(filename,) as f:
        	R_List = pickle.load(f)
	for entry in R_List:
		Histogram.Fill(entry)
	return Histogram
	

def MatchedParticleHistogram(file_path, dR, MomentumThreshold, Histogram=None, BG=False, Save=False, EarlyBreak=0, Flag="", Mode=""):
	print "working on file", file_path
	file = rt.TFile(file_path,"READ")
        # open tree file
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
	if Save and Mode=="decay_vx": R_list = []

        for i in xrange(N):
                if i % 50 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i>=EarlyBreak: break
                tree.GetEntry(i)
		for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
	
			previous_ids = []

                        for k in range(0,tree.nGenParticles):
				if BG:
                                        pdgCriterion = True
                                        statusCriterion = tree.genParticle_status[k] == 2
                                else:
                                        pdgCriterion = abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600
                                        statusCriterion = tree.genParticle_status[k] == 2
                                
                                if statusCriterion and pdgCriterion:
                                        pVector = rt.TLorentzVector()
                                        pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                                                tree.genParticle_phi[k],tree.genParticle_mass[k])
                                        delR = jVector.DeltaR(pVector)
                                        if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
                                                v_p = normalize(np.array([pVector[0], pVector[1], pVector[2]]))
                                                phi = PolarPhi(v_p[0],v_p[1])
                                                theta = Theta(v_p[0],v_p[1],v_p[2])
						
						escape = False		#filter out identical daughters
						if len(previous_ids)>0:
							for prid in previous_ids:
								if abs(abs(prid)-abs(tree.genParticle_pdgId[k])) == 2: escape=True
						if escape: continue
						previous_ids.append(tree.genParticle_pdgId[k])
						R_xy = np.sqrt(tree.genParticle_decayvx_x[k]**2 + tree.genParticle_decayvx_y[k]**2)
                                                if Mode == "decay_vx":
							if Histogram != None: Histogram.Fill(R_xy)
							if Save: R_list.append(R_xy)

	if Save:
		if Mode == "decay_vx":
			filename = "R_List_"+Flag+".pkl"
			with open(filename, 'w') as f:
                        	pickle.dump(R_list, f)
			print "saved as", filename 
	return Histogram

def SearchFor(searchingfor, tree, jVector,dR, MomentumThreshold, BG=False):
	previous_ids = []
	for k in range(0,tree.nGenParticles):
		if BG:
                	pdgCriterion = True
                        statusCriterion = tree.genParticle_status[k] == 2
                else:
	                pdgCriterion = abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600
        		statusCriterion = tree.genParticle_status[k] == 2
			if searchingfor == "quark":
				pdgCriterion = abs(tree.genParticle_pdgId[k])==5 
				statusCriterion = tree.genParticle_status[k] == 23
			if searchingfor == "hadron":
				pdgCriterion = abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600
				statusCriterion = tree.genParticle_status[k] == 2
                        if statusCriterion and pdgCriterion:
	                        pVector = rt.TLorentzVector()
        	                pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                		        tree.genParticle_phi[k],tree.genParticle_mass[k])
                        	delR = jVector.DeltaR(pVector)
                                if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
					escape = False		#filter out identical daughters
					if len(previous_ids)>0:
						for prid in previous_ids:
							if abs(abs(prid)-abs(tree.genParticle_pdgId[k])) == 2: escape=True
					if escape: continue
					previous_ids.append(tree.genParticle_pdgId[k])

					if searchingfor == "quark":
						return (tree.genParticle_pt[k],tree.genParticle_eta[k],tree.genParticle_phi[k],k)
					if searchingfor == "hadron":
						return (tree.genParticle_pt[k],tree.genParticle_eta[k],tree.genParticle_phi[k],k)
	return False

def QuarkHadronComparison(file_path, title, dR, MomentumThreshold, BG=False,save=False, EarlyBreak=0):
	print "working on file", file_path
	file = rt.TFile(file_path,"READ")
        # open tree file
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
	#pt_diff, eta_diff, phi_diff = [],[],[]
	PT_DIFF = rt.TH1D("pt_diff","pt_diff",40,0,700)
        ETA_DIFF = rt.TH1D("eta_diff","eta_diff",40,0,0.2)
        PHI_DIFF = rt.TH1D("phi_diff","phi_diff",40,0,0.2)

        for i in xrange(N):
                if i % 50 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i>=EarlyBreak: break
                tree.GetEntry(i)
		for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
			quark_pt_eta_phi = SearchFor("quark", tree, jVector, dR, MomentumThreshold, BG)
			hadron_pt_eta_phi = SearchFor("hadron", tree, jVector, dR, MomentumThreshold, BG)
			if (quark_pt_eta_phi != False and hadron_pt_eta_phi != False):
				pt_diff_temp = abs(quark_pt_eta_phi[0] - hadron_pt_eta_phi[0])
				eta_diff_temp = abs(quark_pt_eta_phi[1] - hadron_pt_eta_phi[1])
				phi_diff_temp = abs(quark_pt_eta_phi[2] - hadron_pt_eta_phi[2])
				PT_DIFF.Fill(pt_diff_temp)
				ETA_DIFF.Fill(eta_diff_temp)
				PHI_DIFF.Fill(phi_diff_temp)
	
	Histogramize([(PT_DIFF,"pT_diff_"+title)], (0,700), "pt_diff"+title, "pt_diff", "# particles", Save=save,Normalize=False, t_sleep=5)
	Histogramize([(ETA_DIFF,"eta_diff_"+title)], (0,0.2), "eta_diff"+title, "eta_diff", "# particles", Save=save,Normalize=False, t_sleep=5)
	Histogramize([(PHI_DIFF,"phi_diff"+title)], (0,0.2), "phi_diff"+title, "phi_diff", "# particles", Save=save,Normalize=False, t_sleep=5)
	#return (pt_diff,eta_diff,phi_diff)			

def CheckInvariantMass(title, file_paths, EarlyBreak=0):
	M_hist = rt.TH1D("M","M",60,0,5000)
	DeltaPhi_hist = rt.TH1D("DeltaPhi","DeltaPhi",60,0,6.5)

	for file_path in file_paths:
                try:
                        
			print "working on file", file_path
			file = rt.TFile.Open(file_path)
        		# open tree file
        		tree = file.Get("demo/tree")
        		N = tree.GetEntries()
			print "number of events =",N

        		for i in xrange(N):
        		        if i % 100 == 0: print "Working on event " ,i
        		        if EarlyBreak > 0 and i>=EarlyBreak: break
        		        tree.GetEntry(i)
				#print "amount of jets =", tree.nJets
				jetlist = []
				totalJet = rt.TLorentzVector()
				if tree.nJets == 2:
					for j in range(0,tree.nJets):
						jetlist.append(rt.TLorentzVector())
						jetlist[j].SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])	
						totalJet = totalJet + jetlist[j]
		
					M_hist.Fill(totalJet.M())
					if tree.nJets == 2: DeltaPhi_hist.Fill(abs(jetlist[0].DeltaPhi(jetlist[1])))
		except:
			print "skip file because of error"
			continue	

	tfile = rt.TFile("histogram_files/CheckInvariantMass_{}.root".format(title),"recreate")
	M_hist.Write()
        DeltaPhi_hist.Write()

	
	'''
		for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
			quark_pt_eta_phi = SearchFor("quark", tree, jVector, dR, MomentumThreshold, BG)
			hadron_pt_eta_phi = SearchFor("hadron", tree, jVector, dR, MomentumThreshold, BG)
			if (quark_pt_eta_phi != False and hadron_pt_eta_phi != False):
				pt_diff_temp = abs(quark_pt_eta_phi[0] - hadron_pt_eta_phi[0])
				eta_diff_temp = abs(quark_pt_eta_phi[1] - hadron_pt_eta_phi[1])
				phi_diff_temp = abs(quark_pt_eta_phi[2] - hadron_pt_eta_phi[2])
				PT_DIFF.Fill(pt_diff_temp)
				ETA_DIFF.Fill(eta_diff_temp)
				PHI_DIFF.Fill(phi_diff_temp)
	'''

if __name__ == '__main__':

	#Files
	
	Signal_2TeV_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_2000/btagHits_noPU/180628_091836/0000/flatTuple_{}.root' #18-33
        Signal_4TeV_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_noPU/180628_092048/0000/flatTuple_{}.root' #18-33
        Signal_2TeV_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_2000/btagHits_wPU/180604_113337/0000/flatTuple_{}.root' #1-17
        Signal_4TeV_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_wPU/180604_071651/0000/flatTuple_{}.root' #1-16
        BG_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_wPVs/180502_130824/0000/flatTuple_{}.root' #1-47
        BG_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_wPU-v2/180528_101050/0000/flatTuple_{}.root' #1-500
	Signal_2TeV_noPU_String_list=[]
        Signal_4TeV_noPU_String_list=[]
        for i in range(18,34):
                Signal_2TeV_noPU_String_list.append(Signal_2TeV_noPU_String.format(i))
        for i in range(18,34):
                Signal_4TeV_noPU_String_list.append(Signal_4TeV_noPU_String.format(i))

	

	CheckInvariantMass("2TeV",Signal_2TeV_noPU_String_list, EarlyBreak=0)
	CheckInvariantMass("4TeV",Signal_4TeV_noPU_String_list, EarlyBreak=0)


	'''	
	#Deca_vx_r - Histogram
	
	dR = 0.1
	MomentumThreshold = 350
	
	bins = 40
	ran = (0,30)
	title = "Decayvx_R"
	xlabel = "R(x,y)"
	ylabel = "[a.u.]"
	R_Hist_signal1 = rt.TH1D(title, title, bins, ran[0], ran[1])
	R_Hist_signal2 = rt.TH1D(title, title, bins, ran[0], ran[1])
	R_Hist_BG = rt.TH1D(title, title, bins, ran[0], ran[1])
	
	#MatchedParticleHistogram(Signal_file1, dR, MomentumThreshold, Histogram=R_Hist_signal1, BG=False, Save=False, EarlyBreak=100, Flag="Signal1", Mode="decay_vx")
	#MatchedParticleHistogram(Signal_file2, dR, MomentumThreshold, Histogram=R_Hist_signal2, BG=False, Save=False, EarlyBreak=100, Flag="Signal2", Mode="decay_vx")
	#MatchedParticleHistogram(BG_file, dR, MomentumThreshold, Histogram=R_Hist_BG, BG=True, Save=False, EarlyBreak=100, Flag="BG", Mode="decay_vx")
	LoadHistogram(R_Hist_signal1, "decay_vx/R_List_Signal1.pkl")
	LoadHistogram(R_Hist_signal2, "decay_vx/R_List_Signal2.pkl")
	LoadHistogram(R_Hist_BG, "decay_vx/R_List_BG.pkl")
	
	histograms = [(R_Hist_signal1,'2TeV-Signal'), (R_Hist_signal2,'4TeV-Signal'), (R_Hist_BG,'Background')]
	Histogramize(histograms, ran, title, xlabel, ylabel, Save=True,Normalize=True, t_sleep=1000)
	
	
	#pT, eta, phi of b-quark vs B-hadron
	
	#QuarkHadronComparison(Signal_file1,"2TeV-signal", dR, MomentumThreshold, BG=False, save=False, EarlyBreak=1000)






	'''
