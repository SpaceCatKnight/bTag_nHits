import ROOT as rt
from time import sleep
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pickle
import csv
import time
import sys

FeatureDict = {'nEvent':0, 'CSV':1, 'pT_hadron':2, 'pT_jet':3, 'decayvx_R':4, 'L1_0.04':5, 'L2_0.04':6, 'L3_0.04':7, 'L4_0.04':8, 'L1_0.06':9, 'L2_0.06':10, 'L3_0.06':11, 'L4_0.06':12, 'L1_0.08':13, 'L2_0.08':14, 'L3_0.08':15, 'L4_0.08':16, 'L1_0.1':17, 'L2_0.1':18, 'L3_0.1':19, 'L4_0.1':20, 'L1_0.16':21, 'L2_0.16':22, 'L3_0.16':23, 'L4_0.16':24}

#helper functions:

def AngleCorr(angle,border):
        """returns input angle inside period indicated by border"""
        if angle > border:
                return angle-2*border
        elif angle < -border:
                return angle+2*border
        else:
                return angle

def PolarPhi(x,y):
        """returns angle phi according to polar/spherical convention"""
        r = np.sqrt(x**2+y**2)
        if x>0:
                return np.arctan(y/x)
        if x<0 and y >= 0:
                return np.arctan(y/x) + np.pi
        if x<0 and y < 0:
                return np.arctan(y/x) - np.pi
        if x==0 and y > 0:
                return np.pi/2
        if x==0 and y < 0:
                return -np.pi/2

def Theta(x,y,z):
        """returns angle theta according to spherical coordinate convention"""
        return np.pi/2 - np.arctan(z/np.sqrt(x**2+y**2))

def normalize(vector):
        """normalizes input vector to unity"""
        r2 = 0
        for i in vector:
                r2 += i**2
        return vector/np.sqrt(r2)

def Eta(theta):
        """converts the angle theta to the pseudorapidity eta"""
        return -np.log(np.tan(theta/2))

def DeltaR(theta1,theta2,phi1,phi2):
        """returns deltaR according to particle physics convention"""
        deta = Eta(theta1)-Eta(theta2)
        dphi = AngleCorr(phi1-phi2,np.pi)
        return np.sqrt(deta**2 + dphi**2)

def DeltaR_eta(eta1,eta2,phi1,phi2):
        """returns deltaR according to particle physics convention using eta"""
        deta = eta1 - eta2
        dphi = AngleCorr(phi1-phi2,np.pi)
        return np.sqrt(deta**2 + dphi**2)

def ShiftedThetaPhi(x,y,z,dx,dy,dz):
        """returns a tuple containing angles theta and phi with respect to a point (dx,dy,dz)"""
        return (Theta(x-dx,y-dy,z-dz),PolarPhi(x-dx,y-dy))

def TrajectoryLength(theta,v_p):
	"""estimates the optimum trajecory length for a particle propagating at angle theta with velocity modulus v_p"""
        if theta <= 0.079*np.pi or theta >= 0.921*np.pi: #Estimating necessary length for trajectory
                return 60/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)
        elif theta >= 0.224*np.pi and theta <= 0.776*np.pi:
                return 25/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)
        else:
                return 45/np.sqrt(v_p[0]**2+v_p[1]**2+v_p[2]**2)

def Initialize3DPlot(title, xlabel, ylabel, zlabel, grid=None):
        '''initialize 3D-plot'''
        fig = plt.figure(title)
        fig.clf()
        ax = Axes3D(fig)
        ax.clear()
        plt.title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        #plot grid of modules for visual reference
        if grid != None:
                ax.scatter(grid[0],grid[1],grid[2],c='k',s=1,linewidths=0.1)
        return ax

def MakeGrid(file_path, EarlyBreak=0):
	"""takes the path to a root file as input and returns the global cartesian coordinates of all the detUnits hit until event denoted by EarlyBreak and saves them to a pkl-file. Used for visual reference in 3D-plots"""
	#hist = rt.TH1D("det_unit_R","det_unit_R (2D)",80,0,20)
        file = rt.TFile.Open(file_path)
	tree = file.Get("demo/tree")
        N = tree.GetEntries()
        DetUnits = []
        for i in xrange(N):
                if i % 50 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i >= EarlyBreak: break
                tree.GetEntry(i)
		#for k in range(tree.nDetUnits):
		#	hist.Fill(np.sqrt(tree.detUnit_X[k]**2 + tree.detUnit_Y[k]**2))
                for DetUnit in zip(tree.detUnit_X,tree.detUnit_Y,tree.detUnit_Z):
                        if DetUnit not in DetUnits:
                                DetUnits.append(DetUnit)
        X, Y, Z = [], [], []
        for DetUnit in DetUnits:
                X.append(DetUnit[0])
                Y.append(DetUnit[1])
                Z.append(DetUnit[2])
        with open("Grid.pkl", 'w') as f:
                pickle.dump([X,Y,Z], f)
        return X, Y, Z

def PlotTrajectory(vx,p,ax,T_max,res,col,lwidth,lstyle):
	"""plots the trajectory of a particle starting at vertex point vx=(x,y,z) with velocity/momentum vector p=(px,py,pz) on 3D plot ax with length T_max, resolution res, color col, line width lwidth and line style lstyle"""
        Tx,Ty,Tz=[],[],[]
        v_t = normalize(np.array([p[0],p[1],p[2]]))
        for t in np.linspace(0,T_max,res):
                Tx.append(vx[0]+t*v_t[0])
                Ty.append(vx[1]+t*v_t[1])
                Tz.append(vx[2]+t*v_t[2])
        ax.plot(xs=Tx,ys=Ty,zs=Tz,color=col,linewidth=lwidth,linestyle=lstyle)#plots all the tracks

def Efficient_Correlation_Hist2D(title, data, feature_x, feature_y, ran_x, ran_y, bins_x, bins_y, Profiling=False, Save=False):
	"""Creates a 2D histogram of features saved in the cluster matching function (see FeatureDict)"""
	f_x, f_y = FeatureDict[feature_x], FeatureDict[feature_y]
	Hist = rt.TH2I(title,title,bins_x,ran_x[0],ran_x[1],bins_y,ran_y[0],ran_y[1])
	for particle in data:
		Hist.Fill(particle[f_x], particle[f_y])
	canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        Hist.GetXaxis().SetTitle(feature_x)
        Hist.GetYaxis().SetTitle(feature_y)
        Hist.GetYaxis().SetTitleOffset(1.5)
	if Profiling:
		Hist.ProfileX().Draw("COLZ")
	else:
        	Hist.Draw("COLZ")
        if Save: 
		canvas.SaveAs('analytic_plots/correlation_2D_'+title+'.png')
		print 'saved as analytic_plots/correlation_2D_'+title+'.png'
	sleep(10)

def Efficient_Layer_Hist(title,data,dR,minR=0,minPT=0,jet_pT=False,DisplayRatios=True,Save=False):
	"""creates a bar plot of the global amount of matched clusters for each layer and also plots the ratio of consequtive layers. dR should be 0.04, 0.06, 0.08, 0.1 or 0.16. All clusters in the file are taken into account. minR (decay vertex R) and minPT can be used for additional constraint."""
        L1,L2,L3,L4, = 0,0,0,0
	if dR == 0.04:
		dR_tag = 5
	elif dR == 0.06:
		dR_tag = 9
	elif dR == 0.08:
		dR_tag = 13
	elif dR == 0.1:
		dR_tag = 17
	elif dR == 0.16:
		dR_tag = 21
	else:
		print "invalid dR-input"
		return False
	if jet_pT:
		JetOrHadron = 'jet'
		pT_tag = 3
	else:
		JetOrHadron = 'hadron'
		pT_tag = 2
        if minR>0 or minPT>0:
                add = ' with R>'+str(minR)+' and '+JetOrHadron+'_pT>'+str(minPT)
        else:
                add = ''
        for particle in data:
                if particle[4] >= minR and particle[pT_tag] >= minPT:
        		L1 += particle[dR_tag]
			L2 += particle[dR_tag+1]
			L3 += particle[dR_tag+2]
			L4 += particle[dR_tag+3]
	if DisplayRatios:
		fig2, ax2 = plt.subplots(1,2,figsize=(9,5))
        	fig2.suptitle('Hit Clusters per Layer'+add+' inside dR<'+str(dR)+' on '+title+' sample')
        	ax2[0].bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
        	ax2[0].set_ylabel('Clusters')
        	ax2[0].set_xticks([0.5,1.5,2.5,3.5])
        	ax2[0].set_xticklabels(['L1','L2','L3','L4'])
        	ax2[1].bar([0.5,1.5,2.5],[L2/float(L1),L3/float(L2),L4/float(L3)],align='center')
        	ax2[1].set_ylabel('[a.u.]')
        	ax2[1].set_xticks([0.5,1.5,2.5])
        	ax2[1].set_xticklabels(['L2/L1','L3/L2','L4/L3'])
        	plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
        else:
		fig2, ax2 = plt.subplots(1,1,figsize=(5,5))
        	fig2.suptitle(title+add+' in dR<'+str(dR))
        	ax2.bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
        	ax2.set_ylabel('Clusters')
        	ax2.set_xticks([0.5,1.5,2.5,3.5])
        	ax2.set_xticklabels(['L1','L2','L3','L4'])
	if Save:
                if minR>0 or minPT>0:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png'
                else:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+title+'.png'
        plt.show()

def Histogramize(Histograms, ran, title, xlabel, ylabel, Save=False,Normalize=True, t_sleep=0):
	"""Draws multiple histograms neatly on a canvas when given to the function as list of tuples (root histogram, title)."""
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
	"""Puts elements of a list that was previously saved as a pkl file into a given histogram."""
        with open(filename,) as f:
                R_List = pickle.load(f)
        for entry in R_List:
                Histogram.Fill(entry)
        return Histogram

def Efficient_SeparateLayerHist(datalist, ran, dR,minPT=0,jet_pT=True, Save=False):
	"""Creates for each layer and for each dataset given to the datalist as tuple (data, title) a histogram with the amount of clusters in the x-axis. dR should be one of the following values: 0.04, 0.06, 0.08, 0.1, 0.16"""
        if dR == 0.04:
		dR_tag = 5
	elif dR == 0.06:
		dR_tag = 9
	elif dR == 0.08:
		dR_tag = 13
	elif dR == 0.1:
		dR_tag = 17
	elif dR == 0.16:
		dR_tag = 21
	else:
		print "invalid dR-input"
		return False
	if jet_pT:
		JetOrHadron = 'jet'
		pT_tag = 3
	else:
		JetOrHadron = 'hadron'
		pT_tag = 2
	canvas = rt.TCanvas('canvas','canvas',800,800)
        canvas.SetTitle("Matched clusters per layer in dR<"+str(dR))
        rt.gStyle.SetOptStat(0)
        canvas.Divide(2,2,0,0)
        canvas.GetPad(1).SetTitle("Layer 1")
        canvas.GetPad(2).SetTitle("Layer 2")
        canvas.GetPad(3).SetTitle("Layer 3")
        canvas.GetPad(4).SetTitle("Layer 4")
        Hist_list = [[],[],[],[]]
        for n, data in enumerate(datalist):
                Hist_list[0].append(rt.TH1D("L1."+str(n),"L1",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[1].append(rt.TH1D("L2."+str(n),"L2",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[2].append(rt.TH1D("L3."+str(n),"L3",ran[1]-ran[0],ran[0],ran[1]))
                Hist_list[3].append(rt.TH1D("L4."+str(n),"L4",ran[1]-ran[0],ran[0],ran[1]))
		for particle in data[0]:
			if particle[pT_tag] < minPT: continue
                	Hist_list[0][n].Fill(particle[dR_tag])
                	Hist_list[1][n].Fill(particle[dR_tag+1])
                	Hist_list[2][n].Fill(particle[dR_tag+2])
                	Hist_list[3][n].Fill(particle[dR_tag+3])

        for l,layer in enumerate(Hist_list):
                canvas.cd(l+1)
                legend = rt.TLegend(0.9,0.9,0.65,0.75)
                for n,Hist in enumerate(layer):
                        Hist.GetXaxis().SetTitle("# clusters")
                        Hist.GetYaxis().SetTitle('[a.u.]')
                        Hist.GetYaxis().SetTitleOffset(1.5)
                        Hist.SetLineColor(n+2)
                        legend.AddEntry(Hist,datalist[n][1])
                        if n==0:
                                Hist.DrawNormalized()
                        else:
                                Hist.DrawNormalized("SAME")
                legend.Draw()
        if Save: canvas.SaveAs("SeparateLayerHistDR"+str(dR)+".png")
        sleep(10)

def efficient_binned_tagged_jets_hist(datalist,discriminant, discriminant_cuts, CSV_cuts, bins, nbins, Difference=False, mode="pT_hadron",Save=False):
	"""creates a histogram for each dataset given as list of tuples (data, title, range) of all the jets that were b-tagged by passing a given a list of cut values corresponding to the given pT-bins for CSV and a given discriminant versus a feature given as string to 'mode'. The histograms are saved to a root file for further use."""
	title = "binned_tagged_jets_vs_"+mode
	AllJetsHistlist = []
	CSVHistlist = []
	DiscriminantHistlist = []
	if mode == "pT_hadron":
		feature = 2
	elif mode == "pT_jet":
		feature = 3
	elif mode == "decay_vx":
		feature = 4
	for n,data in enumerate(datalist): 
                print "working on",data[1]
		ran = data[2]
		AllJetsHistlist.append(rt.TH1D(data[1]+"_AllJets",data[1]+"_"+title,nbins,ran[0],ran[1]))
		AllJetsHistlist[n].SetLineColor(4)
		CSVHistlist.append(rt.TH1D(data[1]+"_CSV",data[1]+"_"+title,nbins,ran[0],ran[1]))
		CSVHistlist[n].SetLineColor(3)
		DiscriminantHistlist.append(rt.TH1D(data[1]+"_Discriminant",data[1]+"_"+title,nbins,ran[0],ran[1]))
                DiscriminantHistlist[n].SetLineColor(2)
                for particle in data[0]:
			bin_number = bin_selection(particle,bins)
			if bin_number == -100: continue
			AllJetsHistlist[n].Fill(particle[feature])
			if particle[1] >= CSV_cuts[bin_number]: CSVHistlist[n].Fill(particle[feature])
                        if Difference:
                                L = particle[8]-particle[5]
                        else:
                                try:
					L = particle[16]/float(particle[13])
				except ZeroDivisionError:
					continue
			if L >= discriminant_cuts[bin_number]: DiscriminantHistlist[n].Fill(particle[feature])
	canvaslist = []
	legendlist = []
	Tfilelist = []
	for n,data in enumerate(datalist):
		canvaslist.append(rt.TCanvas(data[1]+"_canvas","canvas",600,600))
		canvaslist[n].SetTitle(data[1]+"_"+title)
        	rt.gStyle.SetOptStat(0)
		legendlist.append(rt.TLegend(0.9,0.9,0.65,0.75))
		legendlist[n].AddEntry(AllJetsHistlist[n], "All jets")
		legendlist[n].AddEntry(CSVHistlist[n], "CSV")
		legendlist[n].AddEntry(DiscriminantHistlist[n], discriminant)
        	AllJetsHistlist[n].GetXaxis().SetTitle(mode)
        	AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
        	AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
		AllJetsHistlist[n].Draw()
		CSVHistlist[n].Draw("SAME")
		DiscriminantHistlist[n].Draw("SAME")
		legendlist[n].Draw()
		if Save:
			canvaslist[n].SaveAs(title+"_"+data[1]+discriminant+".png")
			Tfilelist.append(rt.TFile("histogram_files/pT_hists/"+title+"_"+data[1]+discriminant+".root","recreate"))
                	print "saved histogram as histogram_files/pT_hists/"+title+"_"+data[1]+discriminant+".root"
			AllJetsHistlist[n].Write()
			CSVHistlist[n].Write()
			DiscriminantHistlist[n].Write()

def efficient_tagged_jets_hist(datalist,discriminant, discriminant_cut, CSV_cut, bins, Difference=False, mode="pT_hadron",Save=False):
	"""creates a histogram for each dataset given as list of tuples (data, title, range) of all the jets that were b-tagged by passing a given cut value for CSV and a given discriminant versus a feature given as string to 'mode' (see FEatureDict). The histograms are saved to a root file for further use."""
	title = "eff_tagged_jets_vs_"+mode
	AllJetsHistlist = []
	CSVHistlist = []
	DiscriminantHistlist = []
	if mode == "pT_hadron":
		feature = 2
	elif mode == "pT_jet":
		feature = 3
	elif mode == "decay_vx":
		feature = 4
	for n,data in enumerate(datalist): 
                print "working on",data[1]
		ran = data[2]
		AllJetsHistlist.append(rt.TH1D(data[1]+"_AllJets",data[1]+"_"+title,bins,ran[0],ran[1]))
		AllJetsHistlist[n].SetLineColor(4)
		CSVHistlist.append(rt.TH1D(data[1]+"_CSV",data[1]+"_"+title,bins,ran[0],ran[1]))
		CSVHistlist[n].SetLineColor(3)
		DiscriminantHistlist.append(rt.TH1D(data[1]+"_Discriminant",data[1]+"_"+title,bins,ran[0],ran[1]))
                DiscriminantHistlist[n].SetLineColor(2)
                for particle in data[0]:
			AllJetsHistlist[n].Fill(particle[feature])
			if particle[1] >= CSV_cut: CSVHistlist[n].Fill(particle[feature])
                        if Difference:
                                L = particle[8]-particle[5]
                        else:
                                try:
					L = particle[16]/float(particle[13])
				except ZeroDivisionError:
					continue
			if L >= discriminant_cut: DiscriminantHistlist[n].Fill(particle[feature])
	canvaslist = []
	legendlist = []
	Tfilelist = []
	for n,data in enumerate(datalist):
		canvaslist.append(rt.TCanvas(data[1]+"_canvas","canvas",600,600))
		canvaslist[n].SetTitle(data[1]+"_"+title)
        	rt.gStyle.SetOptStat(0)
		legendlist.append(rt.TLegend(0.9,0.9,0.65,0.75))
		legendlist[n].AddEntry(AllJetsHistlist[n], "All jets")
		legendlist[n].AddEntry(CSVHistlist[n], "CSV")
		legendlist[n].AddEntry(DiscriminantHistlist[n], discriminant)
        	AllJetsHistlist[n].GetXaxis().SetTitle(mode)
        	AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
        	AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
		AllJetsHistlist[n].Draw()
		CSVHistlist[n].Draw("SAME")
		DiscriminantHistlist[n].Draw("SAME")
		legendlist[n].Draw()
		if Save:
			canvaslist[n].SaveAs(title+"_"+data[1]+discriminant+".png")
			Tfilelist.append(rt.TFile("histogram_files/pT_hists/"+title+"_"+data[1]+discriminant+".root","recreate"))
                	AllJetsHistlist[n].Write()
			CSVHistlist[n].Write()
			DiscriminantHistlist[n].Write()

def efficient_exclusive_tagged_jets_hist(signal_title, data, discriminant, discriminant_cut, CSV_cut,ran, bins, Difference=False, mode="pT_hadron",Save=False):
	"""creates three histograms when given a signal dataset: -all jets tagged by CSV but not by the discriminant (L4-L1 for Difference = True, L4/L1 for Difference = False inside given dR), -all jets tagged by the discriminant but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (see FeatureDict)"""
	title = signal_title+"_tagged_jets_vs_"+mode+"_exclusive"
	#AllJetsHist = rt.TH1D(signal_title+"_AllJets",title,bins,ran[0],ran[1])
	#AllJetsHist.SetLineColor(4)
	CSV_and_not_Discriminant = rt.TH1D(signal_title+"_CSV_not_Discriminant",title,bins,ran[0],ran[1])
	CSV_and_not_Discriminant.SetLineColor(3)
	Discriminant_and_not_CSV = rt.TH1D(signal_title+"_Discriminant_and_not_CSV",title,bins,ran[0],ran[1])
        Discriminant_and_not_CSV.SetLineColor(2)
	Discriminant_and_CSV = rt.TH1D(signal_title+"_Discriminant_and_CSV",title,bins,ran[0],ran[1])
        Discriminant_and_CSV.SetLineColor(4)
        for particle in data:
		CSV_tag, Disc_tag = False, False
		if particle[1] >= CSV_cut: 
			CSV_tag = True
                if Difference:
                        L = particle[8]-particle[5] 
                else:
			try:
                        	L = particle[8]/particle[5]
			except ZeroDivisionError:
				continue
		if L >= discriminant_cut: 
				Disc_tag = True
		if Disc_tag and not CSV_tag: Discriminant_and_not_CSV.Fill(particle[FeatureDict[mode]]) 
		if CSV_tag and not Disc_tag: CSV_and_not_Discriminant.Fill(particle[FeatureDict[mode]]) 
		if CSV_tag and Disc_tag: Discriminant_and_CSV.Fill(particle[FeatureDict[mode]]) 
	canvas = rt.TCanvas("canvas","canvas",600,600)
	canvas.SetTitle(title)
	rt.gStyle.SetOptStat(0)
	legend = rt.TLegend(0.9,0.9,0.65,0.75)
	legend.AddEntry(Discriminant_and_not_CSV,discriminant+"_and_not_CSV")
	legend.AddEntry(CSV_and_not_Discriminant,"CSV_and_not_"+discriminant)
	legend.AddEntry(Discriminant_and_CSV, discriminant+"_and_CSV")
	Discriminant_and_not_CSV.GetXaxis().SetTitle(mode)
	Discriminant_and_not_CSV.GetYaxis().SetTitle('# jets')
	Discriminant_and_not_CSV.GetYaxis().SetTitleOffset(1.5)
	Discriminant_and_not_CSV.Draw()
	CSV_and_not_Discriminant.Draw("SAME")
	Discriminant_and_CSV.Draw("SAME")
	legend.Draw()
	if Save:
		canvas.SaveAs(title+discriminant+".png")
		Tfile= rt.TFile("histogram_files/pT_hists/"+title+"_"+discriminant+".root","recreate")
        	Discriminant_and_not_CSV.Write()
		CSV_and_not_Discriminant.Write()
		Discriminant_and_CSV.Write()

def RebinHist(hist,name,bins):
	"""Rebins an efficiency-vs-pT histogram such that the lower statistics at high pT are taken into account"""
	import array
	bins_ = array.array('d',bins)
	return 	hist.Rebin(len(bins_)-1,"rebinned_"+name,bins_)

def Efficiency_vs_pT(title,histlist, hist_all_jets,y_max, feature,Save=False,legend_shift=False):
	"""plots for each histogram of tagged jets given in a list of tuples (histogram, title) the efficiency for each bin, where the x-axis corresponds to the feature given as string (see FeatureDict)."""
	canvas = rt.TCanvas('canvas','canvas',600,600)
	if legend_shift:
		legend = rt.TLegend(0.1,0.1,0.35,0.25)
	else:
		legend = rt.TLegend(0.1,0.9,0.35,0.75)
	graphlist = []
	for n,hist in enumerate(histlist):
		graphlist.append(rt.TGraphAsymmErrors())
		if n==0: graphlist[n].SetTitle(title+"_vs_"+feature)
		graphlist[n].Divide(hist[0],hist_all_jets,"cl=0.683 b(1,1) mode")
		legend.AddEntry(graphlist[n], histlist[n][1],"LEP")
		graphlist[n].SetLineColor(n+2)
		if n<1:
	   		graphlist[n].GetXaxis().SetTitle(feature)
        		graphlist[n].GetYaxis().SetTitle('efficiency')
        		graphlist[n].GetYaxis().SetTitleOffset(1.5)
			graphlist[n].SetMinimum(0.)
			graphlist[n].SetMaximum(y_max)
			graphlist[n].Draw()
		else:
			graphlist[n].Draw("SAME")	
	legend.Draw()
	if Save: canvas.SaveAs(title+"_vs_"+feature+".png")

def efficient_Make_ROC_histograms(title, data, pT_cut,Check=False):
	"""uses data made by Efficient_Cluster_Matcher() and creates for L4-L1, L4/L1 and CSV three histograms each: one without pT-Cut, one with the cut on the jet and one with the cut on the hadron. The histograms are saved as root files for reuse by the Make_ROC_Curves() function"""
	diff_ran = (-25,25)
	diff_bins = diff_ran[1]-diff_ran[0]
	ratio_ran = (0,10)
	ratio_bins = 60
	Diff_hist = rt.TH1D("L4-L1","L4-L1",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_hadron_pT = rt.TH1D("L4-L1_pTH","L4-L1_pTH",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_jet_pT = rt.TH1D("L4-L1_pTJ","L4-L1_pTJ",diff_bins,diff_ran[0],diff_ran[1])
	Ratio_hist = rt.TH1D("L4_L1","L4_L1",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_hadron_pT = rt.TH1D("L4_L1_pTH","L4_L1_pTH",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_jet_pT = rt.TH1D("L4_L1_pTJ","L4_L1_pTJ",ratio_bins,ratio_ran[0],ratio_ran[1])
	CSV_hist = rt.TH1D("CSV","CSV",ratio_bins,0,1)
	CSV_hist_hadron_pT = rt.TH1D("CSV_pTH","CSV_pTH",ratio_bins,0,1)
	CSV_hist_jet_pT = rt.TH1D("CSV_pTJ","CSV_pTJ",ratio_bins,0,1)
	
	ZeroDiv, ZeroDiv_hadron, ZeroDiv_jet = 0,0,0
	
	for particle in data:
		Diff_hist.Fill(particle[8]-particle[5])
		CSV_hist.Fill(particle[1])
		try:
			L4_L1 = particle[20]/particle[17]
			Ratio_hist.Fill(L4_L1)
		except ZeroDivisionError:
			ZeroDiv += 1
		if particle[2] >= pT_cut:
			CSV_hist_hadron_pT.Fill(particle[1])
			Diff_hist_hadron_pT.Fill(particle[8]-particle[5])
			try:
				L4_L1 = particle[20]/particle[17]
				Ratio_hist_hadron_pT.Fill(L4_L1)
			except ZeroDivisionError:
				ZeroDiv_hadron += 1
		if particle[3] >= pT_cut:
			CSV_hist_jet_pT.Fill(particle[1])
			Diff_hist_jet_pT.Fill(particle[8]-particle[5])
			try:
				L4_L1 = particle[20]/particle[17]
				Ratio_hist_jet_pT.Fill(L4_L1)
			except ZeroDivisionError:
				ZeroDiv_jet += 1

	tfile = rt.TFile("histogram_files/{}_histograms.root".format(title),"recreate")
	Diff_hist.Write()				
        Diff_hist_hadron_pT.Write()	
        Diff_hist_jet_pT.Write()
        Ratio_hist.Write()
        Ratio_hist_hadron_pT.Write()
	Ratio_hist_jet_pT.Write()
	CSV_hist.Write()
        CSV_hist_hadron_pT.Write()
	CSV_hist_jet_pT.Write()
	csv_file = open("histogram_files/{}_ZeroDiv.csv".format(title),"wb")
	writer = csv.writer(csv_file)
	writer.writerow([ZeroDiv,ZeroDiv_hadron,ZeroDiv_jet])
	csv_file.close()
	print "saved zero division occurences in histogram_files/{}_ZeroDiv.csv".format(title)
	if Check == True:
		Histogramize([(Diff_hist,'NC'),(Diff_hist_hadron_pT,'pTH'+str(pT_cut)),(Diff_hist_jet_pT,'pTJ'+str(pT_cut))], diff_ran, 'L4-L1', 'L4-L1', '#', Save=False,Normalize=False, t_sleep=20)
		Histogramize([(Ratio_hist,'NC'),(Ratio_hist_hadron_pT,'pTH'+str(pT_cut)),(Ratio_hist_jet_pT,'pTJ'+str(pT_cut))], ratio_ran, 'L4_L1', 'L4/L1', '#', Save=False,Normalize=False, t_sleep=20)

def bin_selection(particle,bins):
	"""given the data row of one particle (made by the cluster matching algorithm), it returns in which of the given pT-bins the particle belongs"""
	bin_number = -100
	for n, bin_ in enumerate(bins):
		if n >= len(bins)-1: break
		if particle[3] >= bins[n] and particle[3] < bins[n+1]:
			bin_number =  n
	return bin_number
	
def efficient_Make_Binned_ROC_histograms(title, data, bins):
	"""uses data made by Efficient_Cluster_Matcher() and creates for L4-L1, L4/L1 and CSV one histograms each given (jet-)pT-bin. The histograms are saved as root files for reuse by the Make_ROC_Curves() function"""
	diff_ran = (-25,25)
	diff_bins = diff_ran[1]-diff_ran[0]
	ratio_ran = (0,10)
	ratio_bins = 60

	Diff_hist_list = []
	Ratio_hist_list = []
	CSV_hist_list = []
	ZeroDiv_list = []
	for bin_ in range(len(bins)-1):
		Diff_hist_list.append(rt.TH1D("L4-L1_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"L4-L1_"+str(bins[bin_])+"_"+str(bins[bin_+1]),diff_bins,diff_ran[0],diff_ran[1]))
		Ratio_hist_list.append(rt.TH1D("L4_L1_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"L4_L1_"+str(bins[bin_])+"_"+str(bins[bin_+1]),ratio_bins,ratio_ran[0],ratio_ran[1]))
		CSV_hist_list.append(rt.TH1D("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),ratio_bins,0,1))
		ZeroDiv_list.append(0)
	
	for particle in data:
		bin_number = bin_selection(particle,bins)
		if bin_number == -100: continue
		
		Diff_hist_list[bin_number].Fill(particle[8]-particle[5])
		CSV_hist_list[bin_number].Fill(particle[1])
		try:
			L4_L1 = particle[20]/particle[17]
			Ratio_hist_list[bin_number].Fill(L4_L1)
		except ZeroDivisionError:
			ZeroDiv_list[bin_number] += 1

	tfile = rt.TFile("histogram_files/{}_binned_histograms.root".format(title),"recreate")
	for hist in Diff_hist_list:
		hist.Write()
	for hist in Ratio_hist_list:
                hist.Write()
	for hist in CSV_hist_list:
                hist.Write()
	print "saved histograms in histogram_files/{}_binned_histograms.root".format(title)

	csv_file = open("histogram_files/{}_binned_ZeroDiv.csv".format(title),"wb")
	writer = csv.writer(csv_file)
	writer.writerow(ZeroDiv_list)
	csv_file.close()
	print "saved zero division occurences in histogram_files/{}_binned_ZeroDiv.csv".format(title)

def find_cuts(file_path,bins):
	"""given a pre-made root file containing binned ROC histograms (made by efficient_Make_Binned_ROC_histograms()) and the corresponding list of bins, it returns the 10% mistag rate cut on L4-L1, L4/L1 and CSV for each bin"""
	File = rt.TFile.Open(file_path)
	diff_ran = (-25,25)
        diff_bins = diff_ran[1]-diff_ran[0]
        ratio_ran = (0,10)
        ratio_bins = 60
	for bin_ in range(len(bins)-1):
		print "\n"
		print "cuts for bin {}_{}: \n".format(bins[bin_],bins[bin_+1])
		print "Delta:"
		Get_ROC_Efficiencies(File.Get("L4-L1_"+str(bins[bin_])+"_"+str(bins[bin_+1])),diff_ran,diff_bins,0,print_cut=True)
		print "Ratio:"
		Get_ROC_Efficiencies(File.Get("L4_L1_"+str(bins[bin_])+"_"+str(bins[bin_+1])),ratio_ran,ratio_bins,0,print_cut=True)
		print "CSV:"
		Get_ROC_Efficiencies(File.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),ratio_bins,0,print_cut=True)
		print "\n"

def Make_Binned_ROC_Curves(title,Signal_title,Background_title,bins,log=False):
	#hsv = plt.get_cmap('hsv')
	#color = hsv(np.linspace(0,1.0,len(bins)-1))
	#color = ['b', 'g', 'r', 'c', 'm', 'y']
	color = ['rosybrown','firebrick','red','olivedrab','chartreuse','lightseagreen','deepskyblue','royalblue','navy','blue','darkorchid','mediumvioletred']
	ratio_ran = (0,10)
        ratio_bins = 60	
	
	Signal_ZeroDiv = np.loadtxt("histogram_files/{}_ZeroDiv.csv".format(Signal_title),delimiter=',')
	Signal_file = rt.TFile("histogram_files/{}_histograms.root".format(Signal_title),"READ")
	Background_ZeroDiv = np.loadtxt("histogram_files/{}_ZeroDiv.csv".format(Background_title),delimiter=',')
	Background_file = 		rt.TFile("histogram_files/{}_histograms.root".format(Background_title),"READ")

	#Signal_ratio_eff = np.ndarray((0,ratio_bins+1))
	#Signal_CSV_eff = np.ndarray((0,ratio_bins+1))
	#Background_ratio_eff = np.ndarray((0,ratio_bins+1))
	#Background_CSV_eff = np.ndarray((0,ratio_bins+1))

	plt.figure("ROC")
	plt.clf()

	for bin_ in range(len(bins)-1):
		Ratio_Signal_Eff = Get_ROC_Efficiencies(Signal_file.Get("L4_L1_"+str(bins[bin_])+"_"+str(bins[bin_+1])),ratio_ran,ratio_bins,Signal_ZeroDiv[bin_])
		Ratio_BG_Eff = Get_ROC_Efficiencies(Background_file.Get("L4_L1_"+str(bins[bin_])+"_"+str(bins[bin_+1])),ratio_ran,ratio_bins,Background_ZeroDiv[bin_])
		CSV_Signal_Eff = Get_ROC_Efficiencies(Signal_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),ratio_ran,ratio_bins,0)
		CSV_BG_Eff = Get_ROC_Efficiencies(Background_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),ratio_ran,ratio_bins,0)
		if log:
			plt.semilogy(Ratio_Signal_Eff,Ratio_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"_"+str(bins[bin_+1]))
			plt.semilogy(CSV_Signal_Eff,CSV_BG_Eff, color = color[bin_],linestyle = '-.',)

		else:
			plt.plot(Ratio_Signal_Eff,1-Ratio_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"_"+str(bins[bin_+1]))
			plt.plot(CSV_Signal_Eff,1-CSV_BG_Eff, color = color[bin_],linestyle = '-.',)

	if log:
		plt.semilogy([0,0],[0,0],'k-',label = 'L4/L1')
		plt.semilogy([0,0],[0,0],'k-.',label = 'CSV')
		plt.semilogy([0,1],[0.1,0.1],'k:')
		plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"$\epsilon$_background")
		plt.legend(loc=4)
	else:
		plt.plot([0,0],[0,0],'k-',label = 'L4/L1')
		plt.plot([0,0],[0,0],'k-.',label = 'CSV')
		#plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
		plt.plot([0,1],[0.9,0.9],'k:')
		plt.xlabel(r"$\epsilon$_signal")
		plt.ylabel(r"1-$\epsilon$_background")
		plt.legend(loc=3)
	plt.title(title+"_ROC-Curves")
	
	plt.savefig("ROC/{}_ROC_Curves.png".format(title))
	print "saved as ROC/{}_ROC_Curves.png".format(title)
	plt.show()
	'''
	plt.figure("Log_ROC")
	plt.clf()
	plt.semilogy(Signal_Ratio_eff,Background_Ratio_eff,'b-',label='L4_L1')
	plt.semilogy(Signal_Ratio_eff_jet_pT,Background_Ratio_eff,'b--',label='L4_L1_pTJ'+str(pT_cut))
	plt.semilogy(Signal_CSV_eff,Background_CSV_eff,'g-',label='CSV')
        plt.semilogy(Signal_CSV_eff_jet_pT,Background_CSV_eff,'g--',label='CSV_pTJ'+str(pT_cut))
	plt.semilogy([0,1],[0.1,0.1],'k:',label="10% mistag")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	plt.title("ROC-Curves_log")
	plt.legend(loc=4)
	plt.savefig("ROC/{}_ROC_Curves_log.png".format(title))
	plt.show()
	'''

def Make_ROC_Curves(title,Signal_title,Background_title,diff_ran,ratio_ran,ratio_bins,pT_cut):
	"""Uses the histograms pre-made by efficient_Make_ROC_histograms() to draw ROC-curves in two different representations. The signal and background title should correspond to the ones used in the latter function. The printouts yield which cut corresponds to a 10% mistag rate"""
	diff_bins = diff_ran[1]-diff_ran[0]
	Signal_ZeroDiv = np.loadtxt("histogram_files/{}_ZeroDiv.csv".format(Signal_title),delimiter=',')
	Signal_file = 			rt.TFile("histogram_files/{}_histograms.root".format(Signal_title),"READ")
	Signal_Diff_eff =  		Get_ROC_Efficiencies(Signal_file.Get("L4-L1"),diff_ran,diff_bins,0)
	Signal_Diff_eff_jet_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4-L1_pTJ"),diff_ran,diff_bins,0)
	Signal_Ratio_eff = 		Get_ROC_Efficiencies(Signal_file.Get("L4_L1"),ratio_ran,ratio_bins,Signal_ZeroDiv[0])
	Signal_Ratio_eff_jet_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4_L1_pTJ"),ratio_ran,ratio_bins,Signal_ZeroDiv[2])
	Signal_CSV_eff = 		Get_ROC_Efficiencies(Signal_file.Get("CSV"),(0,1),ratio_bins,0)
	Signal_CSV_eff_jet_pT = 	Get_ROC_Efficiencies(Signal_file.Get("CSV_pTJ"),(0,1),ratio_bins,0)

	Background_ZeroDiv = np.loadtxt("histogram_files/{}_ZeroDiv.csv".format(Background_title),delimiter=',')
	Background_file = 		rt.TFile("histogram_files/{}_histograms.root".format(Background_title),"READ")
	print "L4-L1"
	Background_Diff_eff =  		Get_ROC_Efficiencies(Background_file.Get("L4-L1"),diff_ran,diff_bins,0,print_cut=True)
	print "L4-L1_pTH"
	Background_Diff_eff_hadron_pT = Get_ROC_Efficiencies(Background_file.Get("L4-L1_pTH"),diff_ran,diff_bins,0,print_cut=True)
	print "L4-L1_pTJ"
	Background_Diff_eff_jet_pT = 	Get_ROC_Efficiencies(Background_file.Get("L4-L1_pTJ"),diff_ran,diff_bins,0,print_cut=True)
	print "L4_L1"
	Background_Ratio_eff = 		Get_ROC_Efficiencies(Background_file.Get("L4_L1"),ratio_ran,ratio_bins,Background_ZeroDiv[0],print_cut=True)
	print "L4_L1_pTH"
	Background_Ratio_eff_hadron_pT=	Get_ROC_Efficiencies(Background_file.Get("L4_L1_pTH"),ratio_ran,ratio_bins,Background_ZeroDiv[1],print_cut=True)
	print "L4_L1_pTJ"
	Background_Ratio_eff_jet_pT = 	Get_ROC_Efficiencies(Background_file.Get("L4_L1_pTJ"),ratio_ran,ratio_bins,Background_ZeroDiv[2],print_cut=True)
	print "CSV"
	Background_CSV_eff = 		Get_ROC_Efficiencies(Background_file.Get("CSV"),(0,1),ratio_bins,0,print_cut=True)
        print "CSV_pTH"
	Background_CSV_eff_hadron_pT =	Get_ROC_Efficiencies(Background_file.Get("CSV_pTH"),(0,1),ratio_bins,0,print_cut=True)
	print "CSV_pTJ"
	Background_CSV_eff_jet_pT = 	Get_ROC_Efficiencies(Background_file.Get("CSV_pTJ"),(0,1),ratio_bins,0,print_cut=True)
	plt.figure("ROC")
	plt.clf()
	plt.plot(Signal_Diff_eff,1-Background_Diff_eff,'r-',label='L4-L1')
	plt.plot(Signal_Diff_eff_hadron_pT,1-Background_Diff_eff,'r-.',label='L4-L1_pTH'+str(pT_cut))
	plt.plot(Signal_Diff_eff_jet_pT,1-Background_Diff_eff,'r--',label='L4-L1_pTJ'+str(pT_cut))
	plt.plot(Signal_Ratio_eff,1-Background_Ratio_eff,'b-',label='L4_L1')
	plt.plot(Signal_Ratio_eff_hadron_pT,1-Background_Ratio_eff,'b-.',label='L4_L1_pTH'+str(pT_cut))
	plt.plot(Signal_Ratio_eff_jet_pT,1-Background_Ratio_eff,'b--',label='L4_L1_pTJ'+str(pT_cut))
	plt.plot(Signal_CSV_eff,1-Background_CSV_eff,'g-',label='CSV')
        plt.plot(Signal_CSV_eff_hadron_pT,1-Background_CSV_eff,'g-.',label='CSV_pTH'+str(pT_cut))
        plt.plot(Signal_CSV_eff_jet_pT,1-Background_CSV_eff,'g--',label='CSV_pTJ'+str(pT_cut))
	plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"1-$\epsilon$_background")
	plt.title("ROC-Curves")
	plt.legend(loc=3)
	plt.savefig("ROC/{}_ROC_Curves.png".format(title))
	plt.figure("Log_ROC")
	plt.clf()
	plt.semilogy(Signal_Diff_eff,Background_Diff_eff,'r-',label='L4-L1')
	plt.semilogy(Signal_Diff_eff_hadron_pT,Background_Diff_eff,'r-.',label='L4-L1_pTH'+str(pT_cut))
	plt.semilogy(Signal_Diff_eff_jet_pT,Background_Diff_eff,'r--',label='L4-L1_pTJ'+str(pT_cut))
	plt.semilogy(Signal_Ratio_eff,Background_Ratio_eff,'b-',label='L4_L1')
	plt.semilogy(Signal_Ratio_eff_hadron_pT,Background_Ratio_eff,'b-.',label='L4_L1_pTH'+str(pT_cut))
	plt.semilogy(Signal_Ratio_eff_jet_pT,Background_Ratio_eff,'b--',label='L4_L1_pTJ'+str(pT_cut))
	plt.semilogy(Signal_CSV_eff,Background_CSV_eff,'g-',label='CSV')
        plt.semilogy(Signal_CSV_eff_hadron_pT,Background_CSV_eff,'g-.',label='CSV_pTH'+str(pT_cut))
        plt.semilogy(Signal_CSV_eff_jet_pT,Background_CSV_eff,'g--',label='CSV_pTJ'+str(pT_cut))
	plt.semilogy([0,1],[0.1,0.1],'k:',label="10% mistag")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	plt.title("ROC-Curves_log")
	plt.legend(loc=4)
	plt.savefig("ROC/{}_ROC_Curves_log.png".format(title))
	plt.show()

def Get_ROC_Efficiencies(histogram,ran,nCuts,ZeroDiv,print_cut=False):
	"""Helper function used in Make_ROC_Curves(). Given a discriminant histogram, it finds the cut corresponding most closely to a 10% mistag rate"""
	Cuts = np.linspace(ran[0],ran[1],nCuts+1)
	bin_ran = (histogram.GetXaxis().FindBin(ran[0]),histogram.GetXaxis().FindBin(ran[1]))
	Efficiencies = np.zeros(nCuts+1)
	FullIntegral = histogram.Integral(bin_ran[0],bin_ran[1])
	for n,cut in enumerate(Cuts):
		bin_cut = histogram.GetXaxis().FindBin(cut)
		Efficiencies[n] = histogram.Integral(bin_cut,bin_ran[1])/(FullIntegral+ZeroDiv)
	diff = 1
	closest = 0
	if print_cut:
		for n,eff in enumerate(Efficiencies):
			if abs(eff - 0.1) < diff:
				closest = n
				diff = abs(eff - 0.1)
		print "Mistag rate:",Efficiencies[closest], "corresponding to a cut at", Cuts[closest]
	return Efficiencies
		
def Efficient_Cluster_Matcher(title, file_path, MomentumThreshold, BG=False, EarlyBreak=0, Continue=False, Protocol=False):
	"""Matches hit clusters inside dR-cones (around jet-axis) of the 5 sizes: 0.04, 0.06, 0.08, 0.1, 0.16 to b-jets (BG==False) or non-b-jets (BG==True) of the given root file and saves them as ressource efficient numpy arrays in .npy files.
          
        Inputs:
		title:			title used in the filename
                file_path:              path to root file
                Momentumthreshold:      momentum threshold for particles to be counted
		BG:			if set as True, only non-signal jets are matched
                EarlyBreak:             non zero integer which denotes the number of events after which the algorithm should stop
		Continue:		if set as True, a file with identical parameters will be opened and continued instead of starting from the beginning

        Outputs:
                creates numpy array where each row correspond to a jet with the attributes in the following order: nEvent, CSV-tag, hadron_pT, jet_pT, decay_vx, Li_0.04, Li_0.06, Li_0.08, Li_0.1, Li_0.16, nPU """

	t1 = time.time()
	print "working on file", file_path
	file = rt.TFile.Open(file_path)
        # open tree file
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
	print "There are",N,"events in this file."
	if Protocol:
		csv_file = open("Protocol.csv","wb")
        	writer = csv.writer(csv_file)
		writer.writerow([title,0,N])
	if Continue:
		HitClusters = np.load("matched_clusters/MatchedClusters_{}.npy".format(title))
		N_start = int(HitClusters[-1,0]+1)
		print "successfully loaded matched_clusters/MatchedClusters_{}.npy and continue on event {}".format(title,N_start) 
	else:
        	HitClusters = np.ndarray(shape=(0,26))
        	N_start = 0
	if EarlyBreak != 0:
		if N_start >= EarlyBreak:
			print "EarlyBreak cannot be smaller than the starting event of the file that should be continued!"
			return None
		N = EarlyBreak
	dR_list = np.array([0.04,0.06,0.08,0.1,0.16])
	for i in xrange(N_start,N):
                if i % 100 == 0: print "Working on event " ,i
		if i != 0 and i%10000==0:
			print "saving file - do not abort computation now!" 
			np.save("/afs/cern.ch/user/m/msommerh/CMSSW_9_4_0_patch1/src/bTag_nHits/HitAnalyzer/matched_clusters/MatchedClusters_{}.npy".format(title),HitClusters)	
			print "saved as matched_clusters/MatchedClusters_{}.npy".format(title)
			if Protocol: writer.writerow([title,i,N])
                tree.GetEntry(i)
                for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
			previous_ids = []
                        for k in range(0,tree.nGenParticles):
				if BG:
                                	pdgCriterion = abs(tree.genParticle_pdgId[k]) != 5
                                	statusCriterion = tree.genParticle_status[k]== 23
                                else:
                                        pdgCriterion = (abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600) or (abs(tree.genParticle_pdgId[k]) > 5000 and abs(tree.genParticle_pdgId[k]) < 6000) 
                                        statusCriterion = tree.genParticle_status[k] == 2
                                if statusCriterion and pdgCriterion:
                                        pVector = rt.TLorentzVector()
                                        pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                                                tree.genParticle_phi[k],tree.genParticle_mass[k])
                                        delR = jVector.DeltaR(pVector)
                                        if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
						v_p = normalize(np.array([jVector[0], jVector[1], jVector[2]]))
                                                phi = PolarPhi(v_p[0],v_p[1])
                                                theta = Theta(v_p[0],v_p[1],v_p[2])
							
						escape = False          #filter out identical daughters
                                                if len(previous_ids)>0:
                                                        for prid in previous_ids:
                                                                if (abs(abs(prid[0])-abs(tree.genParticle_pdgId[k])) == 2 and abs(prid[0])>100): 
									escape=True
                                                if escape: continue
						
                                                previous_ids.append((tree.genParticle_pdgId[k],delR))

						JetFeatures = np.array([i,tree.jet_bTag[j],tree.genParticle_pt[k],tree.jet_pt[j],np.sqrt(tree.genParticle_decayvx_x[k]**2+tree.genParticle_decayvx_y[k]**2)])

						HitsPerLayer = np.zeros(20)
                                                
						for nModule,lenModule in enumerate(tree.nClusters): #finding all clusters inside deltaR<dR
                                                        for nCluster in xrange(0,lenModule):
								ClusterTheta,ClusterPhi = ShiftedThetaPhi(tree.cluster_globalx[nModule][nCluster],tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster],tree.PV_x[0],tree.PV_y[0],tree.PV_z[0])
                                                                DR = DeltaR(theta,ClusterTheta,phi,ClusterPhi)
								Layer = tree.detUnit_layer[nModule]
                                                                for ndR,dR in enumerate(dR_list):
									if DR<dR and Layer!=0: HitsPerLayer[4*ndR+Layer-1] += 1
									
						HitClusters = np.vstack((HitClusters,np.concatenate((JetFeatures,HitsPerLayer,[len(tree.PV_x)]))))

        print "Total Number of matched high pt jets:",HitClusters.shape[0]
	print "saving file - do not abort computation now!" 
	np.save("/afs/cern.ch/user/m/msommerh/CMSSW_9_4_0_patch1/src/bTag_nHits/HitAnalyzer/matched_clusters/MatchedClusters_{}.npy".format(title),HitClusters)
	print "Saved as matched_clusters/Matched_Clusters_{}.npy".format(title)
	t2 = time.time()
	print "total processing time: {}s".format(t2-t1)
	if Protocol:
		writer.writerow([title,N,N])
		csv_file.close()
	

def Merge_PU(title):
	Matched_Clusters = np.load("matched_clusters/MatchedClusters_{}.npy".format(title))
	PU = np.load("matched_clusters/nPU_{}.npy".format(title))
	merged_PU = []
	for i in Matched_Clusters[:,0]:
		merged_PU.append(PU[int(i)])
	merged_PU = np.array(merged_PU)
	np.save("matched_clusters/MatchedClusters_{}.npy".format(title),np.hstack((Matched_Clusters,merged_PU.reshape((len(merged_PU),1)))))

def Plot_PU(title):
	hist = rt.TH1I('PU_'+title,'PU_'+title, 80, 0, 80)
	i = 1
	Last_Event = 0
	data = np.ndarray(shape=(0,26))
	while True:
		try:
			a = np.load('MatchedClusters_{}_PU{}.npy'.format(title,i))
			a[:,0] += Last_Event
			data = np.vstack((data,a))
			Last_Event += a.shape[0]
			i += 1
		except:
			break
	PU = data[:,25]
	for event in PU:
		hist.Fill(event)
	'''
	Events = data[:,0]
	m_last = 0
	for n in range(int(Events[-1])):
		for m,event in enumerate(Events[m_last:]):
			if event == n:
				hist.Fill(PU[m])
				m_last = m
				break
	'''
	c1 = rt.TCanvas('canvas', 'canvas', 600,600)
	hist.GetXaxis().SetTitle("# PV")
        hist.GetYaxis().SetTitle('# Events')
        hist.GetYaxis().SetTitleOffset(1.5)
	hist.Draw()
	c1.SaveAs('PU_hist_'+title+'.png')
	'''
	plt.figure()
	plt.plot(Events,PU)
	plt.title('PU in each event')
	plt.xlabel('event')
	plt.ylabel('# PV')
	plt.savefig('nPU_'+title+'.png')
	plt.show()
	'''

def load_data(title,folder,size,old_version=False):
		print 'loading {} files'.format(title)
		if old_version:
			p = 25
		else:
			p = 26
		data = np.ndarray((0,p)) 
		for n in range(1,size+1):
			try:
				temp = np.load(folder+'MatchedClusters_'+title+str(n)+'.npy')
				data = np.vstack((data,temp))
			except:
				print 'file nr {} not found'.format(n)
		print title,'files loaded'
		return data

def efficiency_vs_PU(title, data, Delta_Cut, Ratio_Cut, CSV_Cut, y_max, pT_Cut=200, pT_Mode = "jet"):
	#if bins == 0: bins = ran[1]-ran[0]
	ran = (0,80)
	bins = 80
	if pT_Mode == "jet":
		pT_index = 3
	else:
		pT_index = 2
	
	import array
	middle = range(11,41)
	bins_ = array.array('d',[0.0, 5.0, 7.0, 9.0, 11.0]+range(11,41)+[42.0, 44.0, 46.0, 49.0, 52.0, 55.0, 58.0, 65.0, 80])

	#make histograms of efficiency vs PU
	AllJets_Hist = rt.TH1D("AllJets","AllJets",bins,ran[0],ran[1])
	Delta_Hist = rt.TH1D("Delta","Delta",bins,ran[0],ran[1])
	Ratio_Hist = rt.TH1D("Ratio","Ratio",bins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",bins,ran[0],ran[1])
	for particle in data:
		if particle[pT_index] < pT_Cut: continue
		AllJets_Hist.Fill(particle[25])
		if particle[1] >= CSV_Cut: CSV_Hist.Fill(particle[25])
               	L_D = particle[8]-particle[5]
		if L_D >= Delta_Cut: Delta_Hist.Fill(particle[25])
                try:
			L_R = particle[16]/float(particle[13])
			if L_R >= Ratio_Cut: Ratio_Hist.Fill(particle[25])
		except ZeroDivisionError:
			continue

	AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
	Delta_Hist = Delta_Hist.Rebin(len(bins_)-1,"Delta",bins_)
	Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
	CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)

	#MAke Graphs and draw them
	canvas = rt.TCanvas('canvas','canvas',600,600)
	legend = rt.TLegend(0.1,0.9,0.35,0.75)
	Delta_Graph = rt.TGraphAsymmErrors()
	Ratio_Graph = rt.TGraphAsymmErrors()
	CSV_Graph = rt.TGraphAsymmErrors()
	Delta_Graph.SetTitle(title+"_vs_PU_pT{}{}".format(pT_Mode,pT_Cut))
	Delta_Graph.Divide(Delta_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	Delta_Graph.SetLineColor(2)
	Ratio_Graph.SetLineColor(3)
	CSV_Graph.SetLineColor(4)
	legend.AddEntry(Delta_Graph, "L4-L1", "LEP")
	legend.AddEntry(Ratio_Graph, "L4/L1", "LEP")
	legend.AddEntry(CSV_Graph, "CSV", "LEP")
	Delta_Graph.GetXaxis().SetTitle("PU")
        Delta_Graph.GetYaxis().SetTitle('efficiency')
        Delta_Graph.GetYaxis().SetTitleOffset(1.5)
	Delta_Graph.SetMinimum(0.)
	Delta_Graph.SetMaximum(y_max)
	Delta_Graph.Draw()
	Ratio_Graph.Draw("SAME")	
	CSV_Graph.Draw("SAME")	
	legend.Draw()
	canvas.SaveAs(title+"_vs_PU_pT{}{}.png".format(pT_Mode,pT_Cut))

def binned_efficiency_vs_PU(title, data, Delta_Cuts, Ratio_Cuts, CSV_Cuts, bins, y_max, pT_Cut=200, pT_Mode = "jet"):
	#if nbins == 0: nbins = ran[1]-ran[0]
	ran = (0,80)
	nbins = 80
	if pT_Mode == "jet":
		pT_index = 3
	else:
		pT_index = 2
	
	import array
	bins_ = array.array('d',[0.0, 5.0, 7.0, 9.0, 11.0]+range(11,41)+[42.0, 44.0, 46.0, 49.0, 52.0, 55.0, 58.0, 65.0, 80])

	#make histograms of efficiency vs PU
	AllJets_Hist = rt.TH1D("AllJets","AllJets",nbins,ran[0],ran[1])
	#Delta_Hist = rt.TH1D("Delta","Delta",nbins,ran[0],ran[1])
	Ratio_Hist = rt.TH1D("Ratio","Ratio",nbins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",nbins,ran[0],ran[1])
	for particle in data:
		if particle[pT_index] < pT_Cut: continue
		
		bin_number = bin_selection(particle,bins)
		if bin_number == -100: continue

		AllJets_Hist.Fill(particle[25])
		if particle[1] >= CSV_Cuts[bin_number]: CSV_Hist.Fill(particle[25])
               	#L_D = particle[8]-particle[5]
		#if L_D >= Delta_Cuts[bin_number]: Delta_Hist.Fill(particle[25])
                try:
			L_R = particle[16]/float(particle[13])
			if L_R >= Ratio_Cuts[bin_number]: Ratio_Hist.Fill(particle[25])
		except ZeroDivisionError:
			continue

	AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
	#Delta_Hist = Delta_Hist.Rebin(len(bins_)-1,"Delta",bins_)
	Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
	CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)

	#Make Graphs and draw them
	canvas = rt.TCanvas('canvas','canvas',600,600)
	legend = rt.TLegend(0.1,0.9,0.35,0.75)
	#Delta_Graph = rt.TGraphAsymmErrors()
	Ratio_Graph = rt.TGraphAsymmErrors()
	CSV_Graph = rt.TGraphAsymmErrors()
	Ratio_Graph.SetTitle(title+"_vs_PU_pT{}{}".format(pT_Mode,pT_Cut))
	#Delta_Graph.Divide(Delta_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	#Delta_Graph.SetLineColor(2)
	Ratio_Graph.SetLineColor(3)
	CSV_Graph.SetLineColor(4)
	#legend.AddEntry(Delta_Graph, "L4-L1", "LEP")
	legend.AddEntry(Ratio_Graph, "L4/L1", "LEP")
	legend.AddEntry(CSV_Graph, "CSV", "LEP")
	Ratio_Graph.GetXaxis().SetTitle("#PV")
        Ratio_Graph.GetYaxis().SetTitle('efficiency')
        Ratio_Graph.GetYaxis().SetTitleOffset(1.5)
	Ratio_Graph.SetMinimum(0.)
	Ratio_Graph.SetMaximum(y_max)
	Ratio_Graph.Draw()
	#Delta_Graph.Draw("SAME")	
	CSV_Graph.Draw("SAME")	
	legend.Draw()
	canvas.SaveAs('binned_'+title+"_vs_PU_pT{}{}.png".format(pT_Mode,pT_Cut))



def discriminant_distributions(title, discriminant, Signal_File, BG_File, bins):
	Signal_ratio_hist_list = []
	BG_ratio_hist_list = []
	canvas_list = []
	legend_list = []

	for n in range(len(bins)-1):
		Signal_ratio_hist_list.append(Signal_File.Get(discriminant+'_'+str(bins[n])+'_'+str(bins[n+1])))
		BG_ratio_hist_list.append(BG_File.Get(discriminant+'_'+str(bins[n])+'_'+str(bins[n+1])))
		
		canvas_list.append(rt.TCanvas('canvas'+str(n),'canvas'+str(n),600,600))
		legend_list.append(rt.TLegend(0.78,0.9,0.54,0.79))
		legend_list[n].AddEntry(Signal_ratio_hist_list[n], "Signal")
        	legend_list[n].AddEntry(BG_ratio_hist_list[n], "BG")
		
		#rt.gStyle.SetOptStat(0)
		Signal_ratio_hist_list[n].SetLineColor(2)
		BG_ratio_hist_list[n].SetLineColor(3)

		Signal_ratio_hist_list[n].GetXaxis().SetTitle(discriminant)
        	Signal_ratio_hist_list[n].GetYaxis().SetTitle('[a.u.]')
        	Signal_ratio_hist_list[n].GetYaxis().SetTitleOffset(1.5)

		#Signal_ratio_hist_list[n].GetYaxis().SetRangeUser(0,BG_ratio_hist_list[n].GetMaximum()*1.4)
		#Signal_ratio_hist_list[n].DrawNormalized()
                #BG_ratio_hist_list[n].DrawNormalized('SAME')
                BG_ratio_hist_list[n].DrawNormalized()
		Signal_ratio_hist_list[n].DrawNormalized('SAME')
		legend_list[n].Draw()
		canvas_list[n].SaveAs('discriminant_distributions/bin'+str(bins[n])+'_'+str(bins[n+1])+'_'+discriminant+'_'+title+'.png')

def find_weight_function(Numerator_data, Denominator_data, feature, binsize):

	if feature == "PU":
		feature_nr = 25
		bins = range(0,81)	
	elif feature == "jet_pT":
		feature_nr = 3
		bins = range(0,3001,binsize)
	nbins = len(bins)-1
	import array
	bins_ = array.array('d',bins)
	Numerator_hist = rt.TH1D("Numerator_hist","Numerator_hist",nbins,bins_)
	Denominator_hist = rt.TH1D("Denominator_hist","Denominator_hist",nbins,bins_)
	
	for particle in Numerator_data:
		Numerator_hist.Fill(particle[feature_nr])
	for particle in Denominator_data:
		Denominator_hist.Fill(particle[feature_nr])

	norm = 1
	
	Numerator_hist.Scale(norm/Numerator_hist.Integral())
	Denominator_hist.Scale(norm/Denominator_hist.Integral())
	Ratio_hist = Numerator_hist.Clone()
	Ratio_hist.SetName("Ratio_hist")
	Ratio_hist.SetTitle("Ratio_hist")
	Ratio_hist.Divide(Denominator_hist)
	#Ratio_hist.Divide(BG_hist)
	#Ratio_hist.Scale(norm/Ratio_hist.Integral())
	#Ratio_Graph = rt.TGraphAsymmErrors()
	#Ratio_Graph = rt.TGraph()
	#Ratio_Graph.SetTitle()
	#Ratio_Graph.Divide(BG_hist,Signal_hist,"cl=0.683 b(1,1) mode")
	#Signal_hist.Scale(norm/Signal_hist.Integral())
	#BG_hist.Scale(norm/BG_hist.Integral())

	#tfile = rt.TFile("ratio_hist.root","recreate")
	#Ratio_hist.Write()
	#Signal_hist.Write()
	#BG_hist.Write()

	ratio_dict = {}
	for bin_nr in range(nbins):
		ratio_dict[bins[bin_nr]] = Ratio_hist.GetBinContent(bin_nr+1)

	print ratio_dict
	
	#test
	Denominator_test_hist = rt.TH1D("Denominator_test_hist","Denominator_test_hist",nbins,bins_)
	for particle in Denominator_data:
		if feature == 'PU':
			if particle[feature_nr]>79: continue
		if feature == 'jet_pT':
			if particle[feature_nr]>3000: continue
		Denominator_test_hist.Fill(particle[feature_nr],ratio_dict[(int(particle[feature_nr])/binsize)*binsize])
	
	Denominator_test_hist.Scale(1.1/Denominator_test_hist.Integral())
	
	canvas = rt.TCanvas('canvas','canvas',600,600)
	legend = rt.TLegend(0.78,0.9,0.54,0.79)
	Numerator_hist.SetLineColor(2)
	Denominator_hist.SetLineColor(3)
	Denominator_test_hist.SetLineColor(4)
	legend.AddEntry(Numerator_hist, 'num')
	legend.AddEntry(Denominator_hist, 'denom')
	legend.AddEntry(Denominator_test_hist, 'denom_test')
	
	Denominator_test_hist.Draw()
	Denominator_hist.Draw("SAME")
	Numerator_hist.Draw("SAME")
	#Denominator_test_hist.Draw("SAME")
	legend.Draw()
	canvas.SaveAs('test_hist.png')
	

if __name__ == '__main__':

        MomentumThreshold = 350
	'''
	if len(sys.argv) > 1: infile = sys.argv[1]
	if len(sys.argv) > 2: jobtitle = sys.argv[2]
	
	Efficient_Cluster_Matcher(jobtitle, infile, MomentumThreshold, BG=True, EarlyBreak=0, Continue=False,Protocol=False)
	'''
	'''
	#select file paths	

	SignalFile1 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M2000_GENSIMDIGIRECO_v2.root"
	#SignalFile2 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M4000_GENSIMDIGIRECO_v2.root"

	Additional_Background_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_wPVs/180502_130824/0000/flatTuple_{}.root'
	'''
	#PU_Background_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_wPU-v2/180528_101050/0000/flatTuple_{}.root'	

	
	#pre-process data

	#Efficient_Cluster_Matcher("4TeV-Signal", SignalFile2, MomentumThreshold, BG=False, EarlyBreak=0, Continue=True)
	#Efficient_Cluster_Matcher("Background", Background, MomentumThreshold, BG=True, EarlyBreak=0, Continue=True,Protocol=True)
	#Efficient_Cluster_Matcher("2TeV-Signal", SignalFile1, 350, BG=False, EarlyBreak=0, Continue=False, Protocol=False)


	#need to get GRID-permission first!!!
	#for n in range(0,0):
		#try:
		#	Efficient_Cluster_Matcher("BG_PU"+str(n), PU_Background_String.format(n), MomentumThreshold, BG=True, EarlyBreak=0, Continue=True,Protocol=False)
		#except:
		#	continue
	'''
	for n in range(29,501):
		try:
			Efficient_Cluster_Matcher("BG_PU"+str(n), PU_Background_String.format(n), MomentumThreshold, BG=True, EarlyBreak=0, Continue=False,Protocol=False)
		except:
			continue

	'''


	'''
	#Make_ROC_histograms("4TeV-Signal", Signal2, 0.04, 0.01, 1200,Check=False)
	
	tagged_jets_hist([(Signal1,'2TeV-signal',(350,1200)),(Signal2,'4TeV-signal',(350,2400))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="pT_hadron",Save=True)
	tagged_jets_hist([(Signal1,'2TeV-signal',(0,1500)),(Signal2,'4TeV-signal',(0,2400))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="pT_jet",Save=True)
	tagged_jets_hist([(Signal1,'2TeV-signal',(0,35)),(Signal2,'4TeV-signal',(0,40))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="decay_vx",Save=True)
	
	tagged_jets_hist([(Signal1,'2TeV-signal',(350,1200)),(Signal2,'4TeV-signal',(350,2400))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="pT_hadron",Save=True)
	tagged_jets_hist([(Signal1,'2TeV-signal',(0,1500)),(Signal2,'4TeV-signal',(0,2400))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="pT_jet",Save=True)
	tagged_jets_hist([(Signal1,'2TeV-signal',(0,35)),(Signal2,'4TeV-signal',(0,40))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="decay_vx",Save=True)
	'''
	#tagged_jets_hist([(Background,'Background',(350,2600))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="pT_hadron",Save=True)
	#tagged_jets_hist([(Background,'Background',(0,2600))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="pT_jet",Save=True)
	#tagged_jets_hist([(Background,'Background',(0,0.1))],"L4-L1", 0.04, 5, 0.65, 60, Difference=True, mode="decay_vx",Save=True)
	
	#tagged_jets_hist([(Background,'Background',(350,2600))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="pT_hadron",Save=True)
	#tagged_jets_hist([(Background,'Background',(0,2600))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="pT_jet",Save=True)
	#tagged_jets_hist([(Background,'Background',(0,0.1))],"L4_L1", 0.1, 1.8, 0.65, 60, Difference=False, mode="decay_vx",Save=True)
	'''
	'''
	#exclusive_tagged_jets_hist('4TeV-Signal', Signal2, "L4-L1", 0.04, 5, 0.675, (350,2200), 60, Difference=True, mode="pT_hadron",Save=True)
	'''
	tagged_jets_hist([(Signal2,"4TeV-Signal",(350,2600)),(Background,'Background',(350,2600))],"L4-L1_pT1200", 0.04, 7, 0.716, 60, Difference=True, mode="pT_hadron",Save=True)
	tagged_jets_hist([(Signal2,"4TeV-Signal",(0,2600)),(Background,'Background',(0,2600))],"L4-L1_pT1200", 0.04, 7, 0.716, 60, Difference=True, mode="pT_jet",Save=True)
	tagged_jets_hist([(Signal2,"4TeV-Signal",(350,2600)),(Background,'Background',(350,2600))],"L4_L1_pT1200", 0.1, 2, 0.716, 60, Difference=False, mode="pT_hadron",Save=True)
	tagged_jets_hist([(Signal2,"4TeV-Signal",(0,2600)),(Background,'Background',(0,2600))],"L4_L1_pT1200", 0.1, 2, 0.716, 60, Difference=False, mode="pT_jet",Save=True)
	'''

	'''
	print 'loading files'
	Signal2 = np.load('MatchedClusters_4TeV-Signal.npy')
	Background = np.ndarray((0,25))
	for n in range(1,28):
		if n==14: continue
		Background = np.vstack((Background,np.load('MatchedClusters_BG{}.npy'.format(n))))
	print 'file loaded'
	
	efficient_tagged_jets_hist([(Background, "Background_npy",(350,2600))],"L4-L1", 5, 0.66, 60, Difference=True, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Background, "Background_npy",(0,2600))],"L4-L1", 5, 0.66, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Background, "Background_npy",(350,2600))],"L4_L1", 1.83, 0.65, 60, Difference=False, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Background, "Background_npy",(0,2600))],"L4_L1", 1.83, 0.65, 60, Difference=False, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal2, "4TeV-signal",(350,2600))],"L4-L1", 5, 0.66, 60, Difference=True, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Signal2, "4TeV-signal",(0,2600))],"L4-L1", 5, 0.66, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal2, "4TeV-signal",(350,2600))],"L4_L1", 1.83, 0.65, 60, Difference=False, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Signal2, "4TeV-signal",(0,2600))],"L4_L1", 1.83, 0.65, 60, Difference=False, mode="pT_jet",Save=True)


	#efficient_Make_ROC_histograms("Background", Background, 1200,Check=False)
	#Make_ROC_Curves("Best_Discriminants","4TeV-Signal","Background",(-25,25),(0,10),60,1200)


	'''
	
	



	'''	
	#efficiency vs pT plots

	hadron_file_diff = rt.TFile("histogram_files/pT_hists/tagged_jets_vs_pT_hadron_dR_0.04_4TeV-signalL4-L1.root","READ")
	hadron_AllJetsHist = hadron_file_diff.Get("4TeV-signal_AllJets")
	hadron_CSVHist = hadron_file_diff.Get("4TeV-signal_CSV")
	hadron_DiffHist = hadron_file_diff.Get("4TeV-signal_Discriminant")
	hadron_file_ratio = rt.TFile("histogram_files/pT_hists/tagged_jets_vs_pT_hadron_dR_0.1_4TeV-signalL4_L1.root","READ")
	hadron_RatioHist = hadron_file_ratio.Get("4TeV-signal_Discriminant")

	hadron_AllJetsHist2 = RebinHist(hadron_AllJetsHist,"AllJets")
	hadron_CSVHist2 = RebinHist(hadron_CSVHist,"CSV")
	hadron_DiffHist2 = RebinHist(hadron_DiffHist,"Diff")
	hadron_RatioHist2 = RebinHist(hadron_RatioHist,"Ratio")
	
	Efficiency_vs_pT("4TeV-Signal_Efficiency_rebin",[(hadron_DiffHist2,"L4-L1"),(hadron_RatioHist2,"L4/L1"),(hadron_CSVHist2,"CSV")], hadron_AllJetsHist2,0.6, "pT_hadron",Save=True,legend_shift=True)
	#Efficiency_vs_pT("4TeV-Signal_Efficiency",[(hadron_DiffHist,"L4-L1"),(hadron_RatioHist,"L4/L1"),(hadron_CSVHist,"CSV")], hadron_AllJetsHist,0.5, "pT_hadron",Save=True)

	jet_file_diff = rt.TFile("histogram_files/pT_hists/tagged_jets_vs_pT_jet_dR_0.04_4TeV-signalL4-L1.root","READ")
	jet_AllJetsHist = jet_file_diff.Get("4TeV-signal_AllJets")
	jet_CSVHist = jet_file_diff.Get("4TeV-signal_CSV")
	jet_DiffHist = jet_file_diff.Get("4TeV-signal_Discriminant")
	jet_file_ratio = rt.TFile("histogram_files/pT_hists/tagged_jets_vs_pT_jet_dR_0.1_4TeV-signalL4_L1.root","READ")
	jet_RatioHist = jet_file_ratio.Get("4TeV-signal_Discriminant")

	jet_AllJetsHist2 = RebinHist(jet_AllJetsHist,"AllJets")
	jet_CSVHist2 = RebinHist(jet_CSVHist,"CSV")
	jet_DiffHist2 = RebinHist(jet_DiffHist,"Diff")
	jet_RatioHist2 = RebinHist(jet_RatioHist,"Ratio")

	Efficiency_vs_pT("4TeV-Signal_Efficiency_rebin",[(jet_DiffHist2,"L4-L1"),(jet_RatioHist2,"L4/L1"),(jet_CSVHist2,"CSV")], jet_AllJetsHist2, 0.6, "pT_jet",Save=True,legend_shift=True)
	#Efficiency_vs_pT("4TeV-Signal_Efficiency",[(jet_DiffHist,"L4-L1"),(jet_RatioHist,"L4/L1"),(jet_CSVHist,"CSV")], jet_AllJetsHist, 0.5, "pT_jet",Save=True)
	
	
	
	BG_hadron_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_Background_npyL4-L1.root","READ")
	BG_hadron_AllJetsHist = BG_hadron_file_diff.Get("Background_npy_AllJets")
	BG_hadron_CSVHist = BG_hadron_file_diff.Get("Background_npy_CSV")
	BG_hadron_DiffHist = BG_hadron_file_diff.Get("Background_npy_Discriminant")
	BG_hadron_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_Background_npyL4_L1.root","READ")
	BG_hadron_RatioHist = BG_hadron_file_ratio.Get("Background_npy_Discriminant")
	
	BG_hadron_AllJetsHist2 = RebinHist(BG_hadron_AllJetsHist,"AllJets")
	BG_hadron_CSVHist2 = RebinHist(BG_hadron_CSVHist,"CSV")
	BG_hadron_DiffHist2 = RebinHist(BG_hadron_DiffHist,"Diff")
	BG_hadron_RatioHist2 = RebinHist(BG_hadron_RatioHist,"Ratio")

	Efficiency_vs_pT("Mistag_Rate_large_file",[(BG_hadron_DiffHist,"L4-L1"),(BG_hadron_RatioHist,"L4/L1"),(BG_hadron_CSVHist,"CSV")], BG_hadron_AllJetsHist,0.4, "pT_hadron",Save=True)
	Efficiency_vs_pT("Mistag_Rate_rebin_large_file",[(BG_hadron_DiffHist2,"L4-L1"),(BG_hadron_RatioHist2,"L4/L1"),(BG_hadron_CSVHist2,"CSV")], BG_hadron_AllJetsHist2,0.3, "pT_hadron",Save=True)

	BG_jet_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_Background_npyL4-L1.root","READ")
	BG_jet_AllJetsHist = BG_jet_file_diff.Get("Background_npy_AllJets")
	BG_jet_CSVHist = BG_jet_file_diff.Get("Background_npy_CSV")
	BG_jet_DiffHist = BG_jet_file_diff.Get("Background_npy_Discriminant")
	BG_jet_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_Background_npyL4_L1.root","READ")
	BG_jet_RatioHist = BG_jet_file_ratio.Get("Background_npy_Discriminant")

	BG_jet_AllJetsHist2 = RebinHist(BG_jet_AllJetsHist,"AllJets")
	BG_jet_CSVHist2 = RebinHist(BG_jet_CSVHist,"CSV")
	BG_jet_DiffHist2 = RebinHist(BG_jet_DiffHist,"Diff")
	BG_jet_RatioHist2 = RebinHist(BG_jet_RatioHist,"Ratio")

	Efficiency_vs_pT("Mistag_Rate_large_file",[(BG_jet_DiffHist,"L4-L1"),(BG_jet_RatioHist,"L4/L1"),(BG_jet_CSVHist,"CSV")], BG_jet_AllJetsHist, 0.4, "pT_jet",Save=True)
	Efficiency_vs_pT("Mistag_Rate_rebin_large_file",[(BG_jet_DiffHist2,"L4-L1"),(BG_jet_RatioHist2,"L4/L1"),(BG_jet_CSVHist2,"CSV")], BG_jet_AllJetsHist2, 0.3, "pT_jet",Save=True)

	'''
	
	#PU_study
	
	
	#load data

	#Signal_PU = load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
	#Background_PU = load_data('BG_PU','matched_clusters/BG_PU/',499)
	
	
	'''
	#PU_ROCs
	
	efficient_Make_ROC_histograms("4TeV-Signal_PU", Signal_PU, 1200,Check=False)
	efficient_Make_ROC_histograms("Background_PU", Background_PU, 1200,Check=False)
	
	Make_ROC_Curves("Best_Discriminants_PU","4TeV-Signal_PU","Background_PU",(-25,25),(0,10),60,1200)
	'''
	'''
	#tagged jets hists

	efficient_tagged_jets_hist([(Background_PU, "Background_PU",(350,2600))],"L4-L1", 5, 0.667, 60, Difference=True, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Background_PU, "Background_PU",(0,2600))],"L4-L1", 5, 0.667, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Background_PU, "Background_PU",(350,2600))],"L4_L1", 2.0, 0.667, 60, Difference=False, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Background_PU, "Background_PU",(0,2600))],"L4_L1", 2.0, 0.667, 60, Difference=False, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal_PU, "4TeV-signal_PU",(350,2600))],"L4-L1", 5, 0.667, 60, Difference=True, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Signal_PU, "4TeV-signal_PU",(0,2600))],"L4-L1", 5, 0.667, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal_PU, "4TeV-signal_PU",(350,2600))],"L4_L1", 2.0, 0.667, 60, Difference=False, mode="pT_hadron",Save=True)
	efficient_tagged_jets_hist([(Signal_PU, "4TeV-signal_PU",(0,2600))],"L4_L1", 2.0, 0.667, 60, Difference=False, mode="pT_jet",Save=True)
	
	#efficiency vs pT plots

	hadron_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_4TeV-signal_PUL4-L1.root","READ")
	hadron_AllJetsHist = hadron_file_diff.Get("4TeV-signal_PU_AllJets")
	hadron_CSVHist = hadron_file_diff.Get("4TeV-signal_PU_CSV")
	hadron_DiffHist = hadron_file_diff.Get("4TeV-signal_PU_Discriminant")
	hadron_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_4TeV-signal_PUL4_L1.root","READ")
	hadron_RatioHist = hadron_file_ratio.Get("4TeV-signal_PU_Discriminant")

	hadron_AllJetsHist2 = RebinHist(hadron_AllJetsHist,"AllJets")
	hadron_CSVHist2 = RebinHist(hadron_CSVHist,"CSV")
	hadron_DiffHist2 = RebinHist(hadron_DiffHist,"Diff")
	hadron_RatioHist2 = RebinHist(hadron_RatioHist,"Ratio")
	
	Efficiency_vs_pT("4TeV-Signal_Efficiency_PU",[(hadron_DiffHist2,"L4-L1"),(hadron_RatioHist2,"L4/L1"),(hadron_CSVHist2,"CSV")], hadron_AllJetsHist2,0.6, "pT_hadron",Save=True,legend_shift=True)
	#Efficiency_vs_pT("4TeV-Signal_Efficiency",[(hadron_DiffHist,"L4-L1"),(hadron_RatioHist,"L4/L1"),(hadron_CSVHist,"CSV")], hadron_AllJetsHist,0.5, "pT_hadron",Save=True)

	jet_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_4TeV-signal_PUL4-L1.root","READ")
	jet_AllJetsHist = jet_file_diff.Get("4TeV-signal_PU_AllJets")
	jet_CSVHist = jet_file_diff.Get("4TeV-signal_PU_CSV")
	jet_DiffHist = jet_file_diff.Get("4TeV-signal_PU_Discriminant")
	jet_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_4TeV-signal_PUL4_L1.root","READ")
	jet_RatioHist = jet_file_ratio.Get("4TeV-signal_PU_Discriminant")

	jet_AllJetsHist2 = RebinHist(jet_AllJetsHist,"AllJets")
	jet_CSVHist2 = RebinHist(jet_CSVHist,"CSV")
	jet_DiffHist2 = RebinHist(jet_DiffHist,"Diff")
	jet_RatioHist2 = RebinHist(jet_RatioHist,"Ratio")

	Efficiency_vs_pT("4TeV-Signal_Efficiency_PU",[(jet_DiffHist2,"L4-L1"),(jet_RatioHist2,"L4/L1"),(jet_CSVHist2,"CSV")], jet_AllJetsHist2, 0.6, "pT_jet",Save=True,legend_shift=True)
	#Efficiency_vs_pT("4TeV-Signal_Efficiency",[(jet_DiffHist,"L4-L1"),(jet_RatioHist,"L4/L1"),(jet_CSVHist,"CSV")], jet_AllJetsHist, 0.5, "pT_jet",Save=True)
	
	
	
	BG_hadron_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_Background_PUL4-L1.root","READ")
	BG_hadron_AllJetsHist = BG_hadron_file_diff.Get("Background_PU_AllJets")
	BG_hadron_CSVHist = BG_hadron_file_diff.Get("Background_PU_CSV")
	BG_hadron_DiffHist = BG_hadron_file_diff.Get("Background_PU_Discriminant")
	BG_hadron_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_hadron_Background_PUL4_L1.root","READ")
	BG_hadron_RatioHist = BG_hadron_file_ratio.Get("Background_PU_Discriminant")
	
	BG_hadron_AllJetsHist2 = RebinHist(BG_hadron_AllJetsHist,"AllJets")
	BG_hadron_CSVHist2 = RebinHist(BG_hadron_CSVHist,"CSV")
	BG_hadron_DiffHist2 = RebinHist(BG_hadron_DiffHist,"Diff")
	BG_hadron_RatioHist2 = RebinHist(BG_hadron_RatioHist,"Ratio")

	#Efficiency_vs_pT("Mistag_Rate_large_file",[(BG_hadron_DiffHist,"L4-L1"),(BG_hadron_RatioHist,"L4/L1"),(BG_hadron_CSVHist,"CSV")], BG_hadron_AllJetsHist,0.4, "pT_hadron",Save=True)
	Efficiency_vs_pT("Mistag_Rate_PU",[(BG_hadron_DiffHist2,"L4-L1"),(BG_hadron_RatioHist2,"L4/L1"),(BG_hadron_CSVHist2,"CSV")], BG_hadron_AllJetsHist2,0.3, "pT_hadron",Save=True)

	BG_jet_file_diff = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_Background_PUL4-L1.root","READ")
	BG_jet_AllJetsHist = BG_jet_file_diff.Get("Background_PU_AllJets")
	BG_jet_CSVHist = BG_jet_file_diff.Get("Background_PU_CSV")
	BG_jet_DiffHist = BG_jet_file_diff.Get("Background_PU_Discriminant")
	BG_jet_file_ratio = rt.TFile("histogram_files/pT_hists/eff_tagged_jets_vs_pT_jet_Background_PUL4_L1.root","READ")
	BG_jet_RatioHist = BG_jet_file_ratio.Get("Background_PU_Discriminant")

	BG_jet_AllJetsHist2 = RebinHist(BG_jet_AllJetsHist,"AllJets")
	BG_jet_CSVHist2 = RebinHist(BG_jet_CSVHist,"CSV")
	BG_jet_DiffHist2 = RebinHist(BG_jet_DiffHist,"Diff")
	BG_jet_RatioHist2 = RebinHist(BG_jet_RatioHist,"Ratio")

	#Efficiency_vs_pT("Mistag_Rate_large_file",[(BG_jet_DiffHist,"L4-L1"),(BG_jet_RatioHist,"L4/L1"),(BG_jet_CSVHist,"CSV")], BG_jet_AllJetsHist, 0.4, "pT_jet",Save=True)
	Efficiency_vs_pT("Mistag_Rate_PU",[(BG_jet_DiffHist2,"L4-L1"),(BG_jet_RatioHist2,"L4/L1"),(BG_jet_CSVHist2,"CSV")], BG_jet_AllJetsHist2, 0.3, "pT_jet",Save=True)
	'''
	'''	
	#Efficiency vs PU

	efficiency_vs_PU('4TeV-Signal_Efficiency', Signal_PU, 5, 2.0, 0.667, (0,80), 1.0, bins=0, pT_Cut=0, pT_Mode="jet")
	'''

	#study on pT-dependent cut (binned)
	'''
	Signal_4TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_4TeV-Signal.npy')
	Signal_2TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_2TeV-Signal.npy')
	Signal_both_noPU = np.vstack((Signal_4TeV_noPU,Signal_2TeV_noPU))
	Background_noPU = load_data('BG','matched_clusters/BG_noPU/',31, old_version=True)
	Signal_4TeV_PU = load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
	Signal_2TeV_PU = load_data('2TeV-Signal_PU','matched_clusters/Signal_PU/',19)
	Signal_both_PU = np.vstack((Signal_4TeV_PU,Signal_2TeV_PU))
	Background_PU = load_data('BG_PU','matched_clusters/BG_PU/',499)
	'''
	#bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000] #used for the pT-depending cuts
	#bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 3000] #used for discriminant distributions
	#coarse_bins = [400,600,1200,3000]
	#coarse_bins = [400,600,1200,1800,3000]
	#coarse_bins = [200,1200,3000]
	coarse_bins = [0,1400,1600,1800,3000]

	bins2 = [350.0, 387.5, 425.0, 462.5, 500.0, 537.5, 575.0, 612.5, 650.0, 687.5, 725.0, 762.5, 800.0, 837.5, 875.0, 912.5, 950.0, 987.5, 1025.0, 1100.0, 1200.0, 1400.0, 1600.0, 1800.0, 2100.0, 2500.0, 3000.0] #used for better looking efficiency vs pT plots
	PU_bins = [0.0, 5.0, 7.0, 9.0, 11.0]+range(11,41)+[42.0, 44.0, 46.0, 49.0, 52.0, 55.0, 58.0, 65.0, 80]


	
	delta_cuts_noPU = [3, 5, 6, 6, 7, 8, 8, 9, 9, 10]
	ratio_cuts_noPU = [1.833, 1.833, 2.0, 2.0, 2.167, 2.167, 2.167, 2.167, 2.167, 2.167]
	CSV_cuts_noPU = [0.633, 0.65, 0.667, 0.683, 0.7, 0.733, 0.733, 0.767, 0.783, 0.767]

	delta_cuts_PU = [4,5,6,7,8,8,9,9,10,10]
	ratio_cuts_PU = [2.0,2.0,2.0,2.167,2.167,2.167,2.167,2.167,2.167,2.167]
	CSV_cuts_PU = [0.65,0.667,0.683,0.7,0.717,0.733,0.75,0.783,0.783,0.783]
	
	'''
	#generate weight function for all cases:
	print "noPU_4TeV_withPT"
	find_weight_function(Background_noPU, Signal_4TeV_noPU, 'jet_pT', 5)
	print 'noPU_both_withPT'
	find_weight_function(Background_noPU, Signal_both_noPU, 'jet_pT', 5)
	print 'withPU_4TeV_withPT'
	find_weight_function(Background_PU, Signal_4TeV_PU, 'jet_pT', 5)
	print 'withPU_both_withPT'
	find_weight_function(Background_PU, Signal_both_PU, 'jet_pT', 5)
	print 'withPU_4TeV_withPV'
	find_weight_function(Background_PU, Signal_4TeV_PU, 'PU',1)
	print 'withPU_both_withPV'
	find_weight_function(Background_PU, Signal_both_PU, 'PU',1)
	'''
	
	#make binned ROC curves for all cases
	'''
	efficient_Make_Binned_ROC_histograms('4TeV-Signal_noPU', Signal_4TeV_noPU, coarse_bins)	
	efficient_Make_Binned_ROC_histograms('both-Signal_noPU', Signal_both_noPU, coarse_bins)	
	efficient_Make_Binned_ROC_histograms('4TeV-Signal_PU', Signal_4TeV_PU, coarse_bins)	
	efficient_Make_Binned_ROC_histograms('both-Signal_PU', Signal_both_PU, coarse_bins)	
	
	efficient_Make_Binned_ROC_histograms('BG_noPU', Background_noPU, coarse_bins)	
	efficient_Make_Binned_ROC_histograms('BG_PU', Background_PU, coarse_bins)	
	'''
	Make_Binned_ROC_Curves('4TeV_noPU','4TeV-Signal_noPU_binned','BG_noPU_binned',coarse_bins,log=True)
	Make_Binned_ROC_Curves('both_Signals_noPU','both-Signal_noPU_binned','BG_noPU_binned',coarse_bins,log=True)
	Make_Binned_ROC_Curves('4TeV_PU','4TeV-Signal_PU_binned','BG_PU_binned',coarse_bins,log=True)
	Make_Binned_ROC_Curves('both_Signals_PU','both-Signal_PU_binned','BG_PU_binned',coarse_bins,log=True)

	#no PU	
	'''
	efficient_binned_tagged_jets_hist([(Signal_noPU, "4TeV-Signal_noPU",(0,2600))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, bins, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Signal_noPU, "4TeV-Signal_noPU",(0,2600))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, bins, 60, Difference=False, mode="pT_jet",Save=True)

	efficient_binned_tagged_jets_hist([(Background_noPU, "Background_noPU",(0,2600))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, bins, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Background_noPU, "Background_noPU",(0,2600))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, bins, 60, Difference=False, mode="pT_jet",Save=True)
	'''
	'''
	jet_file_diff = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_4TeV-Signal_noPUL4-L1.root","READ")
	jet_file_ratio = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_4TeV-Signal_noPUL4_L1.root","READ")

	jet_AllJetsHist2 = RebinHist(jet_file_diff.Get("Signal_noPU_AllJets"),"AllJets",bins)
	jet_CSVHist2 = RebinHist(jet_file_diff.Get("Signal_noPU_CSV"),"CSV",bins)
	jet_DiffHist2 = RebinHist(jet_file_diff.Get("Signal_noPU_Discriminant"),"Diff",bins)
	jet_RatioHist2 = RebinHist(jet_file_ratio.Get("Signal_noPU_Discriminant"),"Ratio",bins)

	Efficiency_vs_pT("binned_4TeV-Signal_Efficiency_noPU",[(jet_DiffHist2,"L4-L1"),(jet_RatioHist2,"L4/L1"),(jet_CSVHist2,"CSV")], jet_AllJetsHist2, 0.6, "pT_jet",Save=True,legend_shift=True)
	
	BG_jet_file_diff = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_Background_noPUL4-L1.root","READ")
	BG_jet_file_ratio = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_Background_noPUL4_L1.root","READ")

	BG_jet_AllJetsHist2 = RebinHist(BG_jet_file_diff.Get("Background_noPU_AllJets"),"AllJets",bins)
	BG_jet_CSVHist2 = RebinHist(BG_jet_file_diff.Get("Background_noPU_CSV"),"CSV",bins)
	BG_jet_DiffHist2 = RebinHist(BG_jet_file_diff.Get("Background_noPU_Discriminant"),"Diff",bins)
	BG_jet_RatioHist2 = RebinHist(BG_jet_file_ratio.Get("Background_noPU_Discriminant"),"Ratio",bins)
	
	Efficiency_vs_pT("binned_Mistag_Rate_noPU",[(BG_jet_DiffHist2,"L4-L1"),(BG_jet_RatioHist2,"L4/L1"),(BG_jet_CSVHist2,"CSV")], BG_jet_AllJetsHist2, 0.3, "pT_jet",Save=True,legend_shift=False)
	'''
	
	#with PU
	'''
	efficient_binned_tagged_jets_hist([(Signal_PU, "4TeV-Signal_PU",(0,2600))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, bins, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Signal_PU, "4TeV-Signal_PU",(0,2600))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, bins, 60, Difference=False, mode="pT_jet",Save=True)

	efficient_binned_tagged_jets_hist([(Background_PU, "Background_PU",(0,2600))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, bins, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Background_PU, "Background_PU",(0,2600))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, bins, 60, Difference=False, mode="pT_jet",Save=True)
	'''
	'''
	jet_file_diff = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_4TeV-Signal_PUL4-L1.root","READ")
	jet_file_ratio = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_4TeV-Signal_PUL4_L1.root","READ")

	jet_AllJetsHist2 = RebinHist(jet_file_diff.Get("4TeV-Signal_PU_AllJets"),"AllJets",bins)
	jet_CSVHist2 = RebinHist(jet_file_diff.Get("4TeV-Signal_PU_CSV"),"CSV",bins)
	jet_DiffHist2 = RebinHist(jet_file_diff.Get("4TeV-Signal_PU_Discriminant"),"Diff",bins)
	jet_RatioHist2 = RebinHist(jet_file_ratio.Get("4TeV-Signal_PU_Discriminant"),"Ratio",bins)

	Efficiency_vs_pT("binned_4TeV-Signal_Efficiency_PU",[(jet_DiffHist2,"L4-L1"),(jet_RatioHist2,"L4/L1"),(jet_CSVHist2,"CSV")], jet_AllJetsHist2, 0.6, "pT_jet",Save=True,legend_shift=True)
	
	BG_jet_file_diff = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_Background_PUL4-L1.root","READ")
	BG_jet_file_ratio = rt.TFile("histogram_files/pT_hists/binned_tagged_jets_vs_pT_jet_Background_PUL4_L1.root","READ")

	BG_jet_AllJetsHist2 = RebinHist(BG_jet_file_diff.Get("Background_PU_AllJets"),"AllJets",bins)
	BG_jet_CSVHist2 = RebinHist(BG_jet_file_diff.Get("Background_PU_CSV"),"CSV",bins)
	BG_jet_DiffHist2 = RebinHist(BG_jet_file_diff.Get("Background_PU_Discriminant"),"Diff",bins)
	BG_jet_RatioHist2 = RebinHist(BG_jet_file_ratio.Get("Background_PU_Discriminant"),"Ratio",bins)

	Efficiency_vs_pT("binned_Mistag_Rate_PU",[(BG_jet_DiffHist2,"L4-L1"),(BG_jet_RatioHist2,"L4/L1"),(BG_jet_CSVHist2,"CSV")], BG_jet_AllJetsHist2, 0.3, "pT_jet",Save=True,legend_shift=False)
	

	#efficient_Make_Binned_ROC_histograms('4TeV-Signal_noPU', Signal_noPU, bins)	
	#efficient_Make_Binned_ROC_histograms('4TeV-Signal_PU', Signal_PU, bins)	
	#efficient_Make_Binned_ROC_histograms('BG_noPU', Background_noPU, bins)	
        #efficient_Make_Binned_ROC_histograms('BG_PU', Background_PU, bins)	

	
	BG_noPU_File = rt.TFile.Open('histogram_files/BG_noPU_binned_histograms.root')
	Signal_noPU_File = rt.TFile.Open('histogram_files/4TeV-Signal_noPU_binned_histograms.root')
	BG_PU_File = rt.TFile.Open('histogram_files/BG_PU_binned_histograms.root')
	Signal_PU_File = rt.TFile.Open('histogram_files/4TeV-Signal_PU_binned_histograms.root')

	discriminant_distributions('noPU', 'L4_L1', Signal_noPU_File, BG_noPU_File, bins)
        discriminant_distributions('PU', 'L4_L1', Signal_PU_File, BG_PU_File, bins)
	discriminant_distributions('noPU', 'CSV', Signal_noPU_File, BG_noPU_File, bins)
	discriminant_distributions('PU', 'CSV', Signal_PU_File, BG_PU_File, bins)
	
	
	binned_efficiency_vs_PU('4TeV-Signal_Efficiency', Signal_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, bins, 0.8, pT_Cut=200, pT_Mode="jet")
	binned_efficiency_vs_PU('4TeV-Signal_Efficiency', Signal_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, bins, 0.8, pT_Cut=1000, pT_Mode="jet")

	binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, bins, 0.3, pT_Cut=200, pT_Mode="jet")
	binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, bins, 0.3, pT_Cut=1000, pT_Mode="jet")

	'''




