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
import ClusterMatcher as CM
import TrackMatcher as TM
import FinalClusterMatcher as FCM
import Analysis as Ana

FeatureDict = {'nEvent':0, 'CSV':1, 'pT_hadron':2, 'pT_jet':3, 'decayvx_R':4, 'L1_0.04':5, 'L2_0.04':6, 'L3_0.04':7, 'L4_0.04':8, 'L1_0.06':9, 'L2_0.06':10, 'L3_0.06':11, 'L4_0.06':12, 'L1_0.08':13, 'L2_0.08':14, 'L3_0.08':15, 'L4_0.08':16, 'L1_0.1':17, 'L2_0.1':18, 'L3_0.1':19, 'L4_0.1':20, 'L1_0.16':21, 'L2_0.16':22, 'L3_0.16':23, 'L4_0.16':24}	

def  Sample_Analysis(global_title, file_paths, feature):

	if feature == 'jet_mass':
		title = "jet_mass"
		nbins = 60
		ran = (0,400)
	
	elif feature == 'particle_mass':
		title = 'particle_mass'
		nbins = 60
		ran = (0,200)

	elif feature == 'jet_pT':
		title = 'jet_pT'
		nbins = 60
		ran = (0,3000)
	
	elif feature == 'particle_pT':
		title = 'particle_pT'
		nbins = 60
		ran = (0,3000)

	elif feature == 'decayvx':
		title = "decay_vertex_R"
		nbins = 60
		ran = (0,30)

	elif feature == 'CSV':
		title = "CSV"
		nbins = '60'
		ran = (0,1)	

	else:
		print "unspecified feature"
		return None

	hist = rt.TH1D(title,title,nbins,ran[0],ran[1])
	for n,file_path in enumerate(file_paths):
		try:
			data_file = rt.TFile.Open(file_path)
			tree = data_file.Get("demo/tree")
			N = tree.GetEntries()
			for i in xrange(N):
				if i % 1000 == 0: print "Working on event " ,i
				tree.GetEntry(i)
				if feature == 'jet_mass':
					for k in range(0,tree.nJets):
						hist.Fill(tree.jet_mass[k])	
			        elif feature == 'particle_mass':
					for k in range(0,tree.nGenParticles):
	        	                	hist.Fill(tree.genParticle_mass[k])	
			        elif feature == 'jet_pT':
			   		for k in range(0,tree.nJets):               
                	                	hist.Fill(tree.jet_pt[k])	
				elif feature == 'particle_pT':
					for k in range(0,tree.nGenParticles):
	        	                	hist.Fill(tree.genParticle_pt[k])	
			        elif feature == 'decayvx':
					for k in range(0,tree.nGenparticles):
	        	                	hist.Fill(np.sqrt(tree.genParticle_decayvx_x[k]**2 + tree.genParticle_decayvx_y[k]**2))	
			        elif feature == 'CSV':
					for k in range(0,tree.nJets):
                	                	hist.Fill(tree.jet_bTag[k])	
		except:
			continue
	
	tfile = rt.TFile("Thesis_Plots/root_files/"+global_title+".root","recreate")
	hist.Write()

def  Efficient_Sample_Analysis(global_title, data, ran, feature):

	if feature == 'jet_pT':
		title = 'jet_pT'
		nbins = 60
		ran = (0,2500)
		featurenr = 3
	
	elif feature == 'particle_pT':
		title = 'particle_pT'
		nbins = 60
		ran = (0,3000)
		featurenr = 2

	elif feature == 'decayvx':
		title = "decay_vertex_R"
		nbins = 60
		#ran = (0,30)
		featurenr = 4

	elif feature == 'CSV':
		title = "CSV"
		nbins = 60
		ran = (0,1)
		featurenr = 1	

	else:
		print "unspecified feature"
		return None

	hist = rt.TH1D(title,title,nbins,ran[0],ran[1])
	for entry in data[:,featurenr]:
		hist.Fill(entry)
	
	tfile = rt.TFile("Thesis_Plots/root_files/"+global_title+".root","recreate")
	hist.Write()

def dR_tagger(dR):
	"""given a dR-value, id returns the index of L1 in the corresponding dR-cone in the data provided by the efficient_cluster_matcher() function"""
	if dR == 0.04:
                return 5
        elif dR == 0.06:
                return 9
        elif dR == 0.08:
                return 13
        elif dR == 0.1:
                return 17
        elif dR == 0.16:
                return 21
        else:
                print "invalid dR-input"
                return False


def Efficient_Layer_Hist2(title,data,dR,minR=0,minPT1=200,minPT2=1000,Save=False):
        """creates a bar plot of the global amount of matched clusters for each layer and also plots the ratio of consequtive layers. dR should be 0.04, 0.06, 0.08, 0.1 or 0.16. All clusters in the file are taken into account. minR (decay vertex R) and minPT can be used for additional constraint."""
        L1,L2,L3,L4 = 0,0,0,0
        L1_th,L2_th,L3_th,L4_th = 0,0,0,0
	dR_tag = dR_tagger(dR)
	if dR_tag == False: return False
	for particle in data:
		if particle[4] and particle[3] >= minPT1:
        		L1 += particle[dR_tag]
        		L2 += particle[dR_tag+1]
        		L3 += particle[dR_tag+2]
        		L4 += particle[dR_tag+3]		
	for particle in data:
                if particle[4] >= minR and particle[3] > minPT2:
                        L1_th += particle[dR_tag]
                        L2_th += particle[dR_tag+1]
                        L3_th += particle[dR_tag+2]
                        L4_th += particle[dR_tag+3]

        fig2, ax2 = plt.subplots(2,1,figsize=(4.5,9))
        #fig2.suptitle('Hit Clusters per Layer inside dR<'+str(dR)+' on '+title+' sample')
        ax2[0].bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
	ax2[0].set_title(r"jet $p_T$ > {} GeV".format(minPT1))
        ax2[0].set_ylabel('# clusters')
        ax2[0].set_xticks([0.5,1.5,2.5,3.5])
        ax2[0].set_xticklabels(['L1','L2','L3','L4'])
        ax2[1].bar([0.5,1.5,2.5,3.5],[L1_th, L2_th, L3_th, L4_th],align='center')
	ax2[1].set_title(r"jet $p_T$ > {} GeV".format(minPT2))
        #ax2[1].set_ylabel('[a.u.]')
        ax2[1].set_xticks([0.5,1.5,2.5,3.5])
        ax2[1].set_xticklabels(['L1','L2','L3','L4'])
	ax2[1].set_ylabel('# clusters')
        plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
        if Save:
        	fig2.savefig('Thesis_Plots/HitsPerLayer'+title+'.png')
        	print 'saved as Thesis_Plots/HitsPerLayer'+title+'.png'
        #plt.show()

def General_Make_ROC_histograms(title, data, dR, pT_cut):
	dR_tag = dR_tagger(dR)
	if dR_tag == False: return False

        diff_ran = (-22,22)
        diff_bins = diff_ran[1]-diff_ran[0]
        ratio_ran = (0,10)
        ratio_bins = 60
       
	Ratio_hist_21 = rt.TH1D("L2_L1","L2_L1",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_31 = rt.TH1D("L3_L1","L3_L1",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_41 = rt.TH1D("L4_L1","L4_L1",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_32 = rt.TH1D("L3_L2","L3_L2",ratio_bins,ratio_ran[0],ratio_ran[1])
	Ratio_hist_42 = rt.TH1D("L4_L2","L4_L2",ratio_bins,ratio_ran[0],ratio_ran[1])
       	Ratio_hist_43 = rt.TH1D("L4_L3","L4_L3",ratio_bins,ratio_ran[0],ratio_ran[1])
	
	Diff_hist_21 = rt.TH1D("L2-L1","L2-L1",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_31 = rt.TH1D("L3-L1","L3-L1",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_41 = rt.TH1D("L4-L1","L4-L1",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_32 = rt.TH1D("L3-L2","L3-L2",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_42 = rt.TH1D("L4-L2","L4-L2",diff_bins,diff_ran[0],diff_ran[1])
	Diff_hist_43 = rt.TH1D("L4-L3","L4-L3",diff_bins,diff_ran[0],diff_ran[1])

	CSV_hist = rt.TH1D("CSV","CSV",ratio_bins,0,1)
        
	ZeroDiv_21, ZeroDiv_31, ZeroDiv_41, ZeroDiv_32, ZeroDiv_42, ZeroDiv_43 = 0,0,0,0,0,0 

        for particle in data:
		if particle[3] >= pT_cut:
			CSV_hist.Fill(particle[1])
                	
			Diff_hist_21.Fill(particle[dR_tag+1]-particle[dR_tag])
			Diff_hist_31.Fill(particle[dR_tag+2]-particle[dR_tag])
			Diff_hist_41.Fill(particle[dR_tag+3]-particle[dR_tag])
			Diff_hist_32.Fill(particle[dR_tag+2]-particle[dR_tag+1])
			Diff_hist_42.Fill(particle[dR_tag+3]-particle[dR_tag+1])
			Diff_hist_43.Fill(particle[dR_tag+3]-particle[dR_tag+2])

                	if particle[dR_tag] != 0:
				L2_L1 = particle[dR_tag+1]/particle[dR_tag]
				Ratio_hist_21.Fill(L2_L1)
			else:
				ZeroDiv_21 += 1
			if particle[dR_tag] != 0:	
				L3_L1 = particle[dR_tag+2]/particle[dR_tag]
				Ratio_hist_31.Fill(L3_L1)
			else:
				ZeroDiv_31 += 1
                	if particle[dR_tag] != 0:	   
				L4_L1 = particle[dR_tag+3]/particle[dR_tag]
				Ratio_hist_41.Fill(L4_L1)
			else:
				ZeroDiv_41 += 1
			if particle[dR_tag+1] != 0:		
				L3_L2 = particle[dR_tag+2]/particle[dR_tag+1]
				Ratio_hist_32.Fill(L3_L2)
			else:
				ZeroDiv_32 += 1
			if particle[dR_tag+1] != 0:		
				L4_L2 = particle[dR_tag+3]/particle[dR_tag+1]
				Ratio_hist_42.Fill(L4_L2)
			else:
				ZeroDiv_42 += 1
			if particle[dR_tag+2] != 0:	
				L4_L3 = particle[dR_tag+3]/particle[dR_tag+2]
				Ratio_hist_43.Fill(L4_L3)
			else:
				ZeroDiv_43 += 1
        
	tfile = rt.TFile("Thesis_Plots/root_files/GeneralROCHists_{}_dR{}.root".format(title,dR),"recreate")
	print "created file: Thesis_Plots/root_files/GeneralROCHists_{}_dR{}.root".format(title,dR)
        Ratio_hist_21.Write()
	Ratio_hist_31.Write()
	Ratio_hist_41.Write()
	Ratio_hist_32.Write()
	Ratio_hist_42.Write()
	Ratio_hist_43.Write()
	Diff_hist_21.Write()
        Diff_hist_31.Write()
        Diff_hist_41.Write()
        Diff_hist_32.Write()
        Diff_hist_42.Write()
        Diff_hist_43.Write()
	CSV_hist.Write()
        csv_file = open("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_{}_dR{}.csv".format(title,dR),"wb")
        writer = csv.writer(csv_file)
        writer.writerow([ZeroDiv_21, ZeroDiv_31, ZeroDiv_41, ZeroDiv_32, ZeroDiv_42, ZeroDiv_43])
        csv_file.close()
        print "saved zero division occurences in Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_{}_dR{}.csv".format(title,dR)

def General_Make_ROC_Curves(title, histlist,log=False, print_cut=False, dR=None):
	'''histlist entry: (signal_hist, bg_hist, title, diff( True or False), signal_ZeroDiv, bg_ZeroDiv,linestyle(optional))'''
        #hsv = plt.get_cmap('hsv')
        #color = hsv(np.linspace(0,1.0,len(bins)-1))
        #color = ['b', 'g', 'r', 'c', 'm', 'y']
        if len(histlist)<=6:
                color = ['red','green','blue','orange','brown','black']
	#if len(histlist)<=7:
        #        color = ['red','green','blue','orange','brown','cyan','black']

        else:
                color = ['deepskyblue','rosybrown','olivedrab','royalblue','firebrick','chartreuse','navy','red','darkorchid','lightseagreen','mediumvioletred','blue']
        ratio_ran = (0,10)
        ratio_bins = 60 
 	diff_ran = (-22,22)      
        diff_bins = diff_ran[1]-diff_ran[0] 
 
        plt.figure("ROC")
        plt.clf()

        for n,hist in enumerate(histlist):
		if print_cut: print hist[2]
		if hist[3] == "diff":
			Signal_Eff = FCM.Get_ROC_Efficiencies(hist[0],diff_ran,diff_bins,0)
                	BG_Eff = FCM.Get_ROC_Efficiencies(hist[1],diff_ran,diff_bins,0,print_cut=print_cut)
		elif hist[3] == "ratio":
                	Signal_Eff = FCM.Get_ROC_Efficiencies(hist[0],ratio_ran,ratio_bins,hist[4])
                	BG_Eff = FCM.Get_ROC_Efficiencies(hist[1],ratio_ran,ratio_bins,hist[5],print_cut=print_cut)
		elif hist[3] == "CSV":
                	Signal_Eff = FCM.Get_ROC_Efficiencies(hist[0],(0,1),ratio_bins,hist[4])
                	BG_Eff = FCM.Get_ROC_Efficiencies(hist[1],(0,1),ratio_bins,hist[5],print_cut=print_cut)
		else:
			print "unspecified discriminant type in histogram nr {}".format(n)
			return False
                if log:
                	if len(hist) <= 6:
				plt.semilogy(Signal_Eff,BG_Eff, color = color[n], linestyle = '-',label=hist[2])
                	elif len(hist) >= 8:
				plt.semilogy(Signal_Eff,BG_Eff, color = hist[6], linestyle = hist[7],label=hist[2])
			else:
				print "size of histlist entries corrupted"
				return False
		else:
			if len(hist) <= 6:
                        	plt.plot(Signal_Eff,1-BG_Eff, color = color[n], linestyle = '-',label=hist[2])
			elif len(hist) >= 8:
				plt.plot(Signal_Eff,1-BG_Eff, color = hist[6], linestyle = hist[7],label=hist[2])
			else:                                             	
                        	print "size of histlist entries corrupted"
                        	return False
        if log:
                #plt.semilogy([0,0],[0,0],'k-',label = 'L4/L1')
                #plt.semilogy([0,0],[0,0],'k-.',label = 'CSV')
                plt.semilogy([0,1],[0.1,0.1],'k:')
                plt.xlabel(r"signal efficiency")
                plt.ylabel(r"mistag rate")
		plt.ylim(10**(-3),1)
		if dR != None: plt.figtext(0.16,0.83,r'$\Delta R$<{}'.format(dR))
                plt.legend(loc=4)
        else:
                #plt.plot([0,0],[0,0],'k-',label = 'L4/L1')
                #plt.plot([0,0],[0,0],'k-.',label = 'CSV')
                #plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
                plt.plot([0,1],[0.9,0.9],'k:')
                plt.xlabel(r"signal efficiency")
                plt.ylabel(r"purity")
		if dR != None: plt.figtext(0.16,0.83,r'$\Delta R$<{}'.format(dR))
                plt.legend(loc=3)
        #plt.title(title+"_ROC-Curves")

        plt.savefig("Thesis_Plots/{}_ROC_Curves.png".format(title))
        print "saved as Thesis_Plots/{}_ROC_Curves.png".format(title)

def efficient_tagged_jets_hist(datalist,discriminant, discriminant_cut, CSV_cut, bins, Difference=False, mode="pT_jet",Save=False):
        """creates a histogram for each dataset given as list of tuples (data, title, range) of all the jets that were b-tagged by passing a given cut value for CSV and a given discriminant versus a feature given as string to 'mode' (see FeatureDict). The histograms are saved to a root file for further use."""
        title = "tagged_jets_vs_"+mode
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
                                if particle[13] != 0:
                                        L = particle[16]/float(particle[13])
                                else:
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
                #AllJetsHistlist[n].Draw()
                #CSVHistlist[n].Draw("SAME")
                #DiscriminantHistlist[n].Draw("SAME")
                #legendlist[n].Draw()
                if Save:
                        #canvaslist[n].SaveAs("Thesis_Plots/"+title+"_"+data[1]+discriminant+".png")
                        Tfilelist.append(rt.TFile("Thesis_Plots/root_files/"+title+"_"+data[1]+discriminant+".root","recreate"))
                        AllJetsHistlist[n].Write()
                        CSVHistlist[n].Write()
                        DiscriminantHistlist[n].Write()

def efficiency_vs_cut(title, signal_hist, bg_hist, signal_ZeroDiv, bg_ZeroDiv, ran, nCuts):
	Cuts = np.linspace(ran[0],ran[1],nCuts+1)	
	signal_eff = FCM.Get_ROC_Efficiencies(signal_hist,ran,nCuts,signal_ZeroDiv)
	bg_eff = FCM.Get_ROC_Efficiencies(bg_hist,ran,nCuts,bg_ZeroDiv)
	plt.figure()
        plt.clf()
	plt.plot(Cuts,signal_eff,'.',label=r'signal efficiency')
	plt.plot(Cuts,1-bg_eff,'.',label=r'purity')
        if title == "L4_L1":
		plt.xlabel("L4/L1")
	elif title == "L4-L1":
		plt.xlabel("L4-L1")
	elif title == "CSV":
		plt.xlabel("p(b-jet)")
	else:
		print "unidentified discriminant"
		plt.xlabel("cut value")
        plt.ylabel(r"$\epsilon$")
        plt.legend(loc=5)
        plt.savefig("Thesis_Plots/efficiency_vs_cut_{}.png".format(title))
        print "saved as Thesis_Plots/efficiency_vs_cut_{}.png".format(title)

def Efficiency_vs_pT(title,histlist, hist_all_jets,y_max,Save=False,legend_shift=False,BG=False, LargeLegend=False):
        """plots for each histogram of tagged jets given in a list of tuples (histogram, title, colorindex(optional)) the efficiency for each bin, where the x-axis corresponds to the feature given as string (see FeatureDict)."""
        canvas = rt.TCanvas('canvas','canvas',600,600)
        if legend_shift:
		if LargeLegend:
			#legend = rt.TLegend(0.1,0.1,0.4,0.3)
			legend = rt.TLegend(0.6,0.9,0.7,0.9)
		else:
                	#legend = rt.TLegend(0.1,0.1,0.35,0.25)
			legend = rt.TLegend(0.65,0.75,0.9,0.9)
        else:
		if LargeLegend:
			legend = rt.TLegend(0.1,0.9,0.4,0.7)
		else:
                	legend = rt.TLegend(0.1,0.9,0.35,0.75)
        graphlist = []
        for n,hist in enumerate(histlist):
                graphlist.append(rt.TGraphAsymmErrors())
                #if n==0: graphlist[n].SetTitle(title+"_vs_jet-pT")
                graphlist[n].Divide(hist[0],hist_all_jets,"cl=0.683 b(1,1) mode")
                legend.AddEntry(graphlist[n], histlist[n][1],"LEP")
		if len(hist) > 2:
			graphlist[n].SetLineColor(hist[2])
		else:
			if n < 3:
        	        	graphlist[n].SetLineColor(n+2)
			else:
				graphlist[n].SetLineColor(n+3)
                if n<1:
                        graphlist[n].GetXaxis().SetTitle("jet p_{T} (GeV)")
                        if BG:
				graphlist[n].GetYaxis().SetTitle('mistag rate')
			else:
				graphlist[n].GetYaxis().SetTitle('efficiency')
                        graphlist[n].GetYaxis().SetTitleOffset(1.5)
                        graphlist[n].SetMinimum(0.)
                        graphlist[n].SetMaximum(y_max)
                        graphlist[n].Draw()
                else:
                        graphlist[n].Draw("SAME")       
        legend.Draw()
        if Save: canvas.SaveAs("Thesis_Plots/"+title+"_vs_jet-pT.png")

def LoadFromRootFile(file_path, object_name_list):
	tfile = rt.TFile.Open(file_path)
	if len(object_name_list) == 1:
		return tfile.Get(object_name_list[0])
	else:
		objects = []
		for obj in object_name_list:
			objects.append(tfile.Get(obj))
		return objects

def efficient_Make_Binned_ROC_histograms(title, data, bins, PU_range='full'):
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
                if PU_range != 'full':
                        if particle[-1]<PU_range[0] or particle[-1]>PU_range[1]: continue
                bin_number = FCM.bin_selection(particle,bins)
                if bin_number == -100: continue

                Diff_hist_list[bin_number].Fill(particle[8]-particle[5])
                CSV_hist_list[bin_number].Fill(particle[1])
                if particle[17] != 0:
                        L4_L1 = particle[20]/particle[17]
                        Ratio_hist_list[bin_number].Fill(L4_L1)
                else:
                        ZeroDiv_list[bin_number] += 1

        tfile = rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(title),"recreate")
        for hist in Diff_hist_list:
                hist.Write()
        for hist in Ratio_hist_list:
                hist.Write()
        for hist in CSV_hist_list:
                hist.Write()
        print "saved histograms in Thesis_Plots/root_files/{}_histograms.root".format(title)

        csv_file = open("Thesis_Plots/root_files/{}_ZeroDiv.csv".format(title),"wb")
        writer = csv.writer(csv_file)
        writer.writerow(ZeroDiv_list)
        csv_file.close()
        print "saved zero division occurences in Thesis_Plots/root_files/{}_ZeroDiv.csv".format(title)

def Make_Binned_ROC_Curves(title,Signal_title,Background_title,bins, diff=False,log=False,dR=None, ANN=False):
        #hsv = plt.get_cmap('hsv')
        #color = hsv(np.linspace(0,1.0,len(bins)-1))
        #color = ['b', 'g', 'r', 'c', 'm', 'y']
        if len(bins)<=6:
                color = ['red','green','blue','orange','brown']
        else:
                color = ['deepskyblue','rosybrown','olivedrab','royalblue','firebrick','chartreuse','navy','red','darkorchid','lightseagreen','mediumvioletred','blue']
        ratio_ran = (0,10)
        ratio_bins = 60
	diff_ran = (-22,22)
	diff_bins = diff_ran[1]-diff_ran[0]

	if diff:
		nbins = diff_bins
		ran = diff_ran
		dis_string = "L4-L1_"
	else:
		nbins = ratio_bins
		ran = ratio_ran
		dis_string = "L4_L1_"
	
	if ANN:
		nbins = 60
		ran = (0,1)
		dis_string = "ANN_"
	
	if ANN:
		Signal_ZeroDiv = np.zeros(len(bins)-1)
		Background_ZeroDiv = np.zeros(len(bins)-1)
	else:
	        Signal_ZeroDiv = np.loadtxt("Thesis_Plots/root_files/{}_ZeroDiv.csv".format(Signal_title),delimiter=',')
		Background_ZeroDiv = np.loadtxt("Thesis_Plots/root_files/{}_ZeroDiv.csv".format(Background_title),delimiter=',')
        Signal_file = rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(Signal_title),"READ")
        Background_file =    rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(Background_title),"READ")

        plt.figure("ROC")
        plt.clf()

        for bin_ in range(len(bins)-1):
                Dis_Signal_Eff = FCM.Get_ROC_Efficiencies(Signal_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),ran,nbins,Signal_ZeroDiv[bin_])
                Dis_BG_Eff = FCM.Get_ROC_Efficiencies(Background_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),ran,nbins,Background_ZeroDiv[bin_])
                CSV_Signal_Eff = FCM.Get_ROC_Efficiencies(Signal_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),ratio_bins,0)
                CSV_BG_Eff = FCM.Get_ROC_Efficiencies(Background_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),ratio_bins,0)
		if log:
			if bins[bin_] == 0:
				plt.semilogy(Dis_Signal_Eff,Dis_BG_Eff, color = color[bin_], linestyle = '-',label=str(200)+"-"+str(bins[bin_+1])+' GeV')
			else:
                        	plt.semilogy(Dis_Signal_Eff,Dis_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"-"+str(bins[bin_+1])+' GeV')
                        plt.semilogy(CSV_Signal_Eff,CSV_BG_Eff, color = color[bin_],linestyle = '--',)

                else:
			if bins[bin_] == 0:
                        	plt.plot(Dis_Signal_Eff,1-Dis_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"-"+str(bins[bin_+1])+" GeV")
                        else:
				plt.plot(CSV_Signal_Eff,1-CSV_BG_Eff, color = color[bin_],linestyle = '--',)

        if log:
		if ANN:
			plt.semilogy([0,0],[0,0],'k-',label = 'ANN')
		else:
			if diff:
				plt.semilogy([0,0],[0,0],'k-',label = 'L4-L1')
			else:
                		plt.semilogy([0,0],[0,0],'k-',label = 'L4/L1')
                plt.semilogy([0,0],[0,0],'k-.',label = 'CSV')
                plt.semilogy([0,1],[0.1,0.1],'k:')
                plt.xlabel(r"signal efficiency")
                plt.ylabel(r"mistag rate")
		plt.ylim(10**(-3),1)
		if dR != None: plt.figtext(0.16,0.83,r'$\Delta R$<{}'.format(dR))
                plt.legend(loc=4)
        else:
		if ANN:
			plt.plot([0,0],[0,0],'k-',label = 'ANN')
		else:
			if diff:
				plt.plot([0,0],[0,0],'k-',label = 'L4-L1')
			else:
                		plt.plot([0,0],[0,0],'k-',label = 'L4/L1')
                plt.plot([0,0],[0,0],'k-.',label = 'CSV')
                #plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
                plt.plot([0,1],[0.9,0.9],'k:')
                plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"1-$\epsilon$_background")
		if dR != None: plt.figtext(0.16,0.83,r'$\Delta R$<{}'.format(dR))
                plt.legend(loc=3)
        #plt.title(title+"_ROC-Curves")

        plt.savefig("Thesis_Plots/{}_ROC_Curves.png".format(title))
        print "saved as Thesis_Plots/{}_ROC_Curves.png".format(title)

def efficient_binned_tagged_jets_hist(datalist,discriminant, discriminant_cuts, CSV_cuts, bins, nbins, Difference=False, mode="pT_jet",Save=False):
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
                        bin_number = FCM.bin_selection(particle,bins)
                        if bin_number == -100: continue
                        AllJetsHistlist[n].Fill(particle[feature])
                        if particle[1] >= CSV_cuts[bin_number]: CSVHistlist[n].Fill(particle[feature])
                        if Difference:
                                L = particle[8]-particle[5]
                        else:
                                if particle[17] != 0:
                                        L = particle[20]/float(particle[17])
                                else:
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
                AllJetsHistlist[n].GetXaxis().SetTitle("jet p_{T} (GeV)")
                AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
                AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
                #AllJetsHistlist[n].Draw()
                #CSVHistlist[n].Draw("SAME")
                #DiscriminantHistlist[n].Draw("SAME")
                #legendlist[n].Draw()
                if Save:
                        #canvaslist[n].SaveAs(title+"_"+data[1]+discriminant+".png")
                        Tfilelist.append(rt.TFile("Thesis_Plots/root_files/"+title+"_"+data[1]+discriminant+".root","recreate"))
                        print "saved histogram as Thesis_Plots/root_files/"+title+"_"+data[1]+discriminant+".root"
                        AllJetsHistlist[n].Write()
                        CSVHistlist[n].Write()
                        DiscriminantHistlist[n].Write()

def PU_histograms(Signal_2TeV, Signal_4TeV, BG):
	hist_2TeV = rt.TH1D("PU_2TeV","PU_2TeV",80,0,80)
	hist_4TeV = rt.TH1D("PU_4TeV","PU_4TeV",80,0,80)
	hist_BG = rt.TH1D("PU_BG","PU_BG",80,0,80)

	for event in Signal_2TeV[:,25]:
		hist_2TeV.Fill(event)
	for event in Signal_4TeV[:,25]:
		hist_4TeV.Fill(event)
	for event in BG[:,25]:
		hist_BG.Fill(event)
	
	tfile = rt.TFile("Thesis_Plots/root_files/PU_distributions.root","recreate")
	hist_2TeV.Write()
	hist_4TeV.Write()
	hist_BG.Write()

def binned_efficiency_vs_PU(title, data, Delta_Cuts, Ratio_Cuts, CSV_Cuts, bins, y_max, pT_Cut=200, pT_Mode = "jet",BG=False):
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
        Delta_Hist = rt.TH1D("Delta","Delta",nbins,ran[0],ran[1])
        Ratio_Hist = rt.TH1D("Ratio","Ratio",nbins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",nbins,ran[0],ran[1])
        for particle in data:
                if particle[pT_index] < pT_Cut: continue

                bin_number = FCM.bin_selection(particle,bins)
                if bin_number == -100: continue

                AllJets_Hist.Fill(particle[25])
                if particle[1] >= CSV_Cuts[bin_number]: CSV_Hist.Fill(particle[25])
                L_D = particle[8]-particle[5]
                if L_D >= Delta_Cuts[bin_number]: Delta_Hist.Fill(particle[25])
                if particle[17] != 0:
                        L_R = particle[20]/float(particle[17])
                        if L_R >= Ratio_Cuts[bin_number]: Ratio_Hist.Fill(particle[25])
                else:
                        continue

        AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        Delta_Hist = Delta_Hist.Rebin(len(bins_)-1,"Delta",bins_)
        Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)

        #Make Graphs and draw them
        canvas = rt.TCanvas('canvas','canvas',600,600)
	rt.gStyle.SetOptTitle(0)
        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        Delta_Graph = rt.TGraphAsymmErrors()
        Ratio_Graph = rt.TGraphAsymmErrors()
        CSV_Graph = rt.TGraphAsymmErrors()
        #Ratio_Graph.SetTitle(title+"_vs_PU_pT{}{}".format(pT_Mode,pT_Cut))
        Delta_Graph.Divide(Delta_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Delta_Graph.SetLineColor(3)
        Ratio_Graph.SetLineColor(2)
        CSV_Graph.SetLineColor(4)
        legend.AddEntry(Ratio_Graph, "L4/L1", "LEP")
	legend.AddEntry(Delta_Graph, "L4-L1", "LEP")
        legend.AddEntry(CSV_Graph, "CSV", "LEP")
        Ratio_Graph.GetXaxis().SetTitle("#PV")
        if BG:
		Ratio_Graph.GetYaxis().SetTitle('mistag rate')
	else:	
		Ratio_Graph.GetYaxis().SetTitle('efficiency')
        Ratio_Graph.GetYaxis().SetTitleOffset(1.5)
	Ratio_Graph.SetMinimum(0.)
        Ratio_Graph.SetMaximum(y_max)
        Ratio_Graph.Draw()
        Delta_Graph.Draw("SAME")
        CSV_Graph.Draw("SAME")
        legend.Draw()
        canvas.SaveAs('Thesis_Plots/'+title+"_vs_PU_pT{}{}.png".format(pT_Mode,pT_Cut))

def QuarkHadronComparison(file_paths, title, dR, MomentumThreshold, BG=False, EarlyBreak=0):
	"""Makes a histogram of the fractional difference of pT, eta and phi between hadrons and quarks in the sample given. """      
	PT_DIFF = rt.TH1D("pt_diff","pt_diff",40,0,1)
       	ETA_DIFF = rt.TH1D("eta_diff","eta_diff",40,0,1)
       	PHI_DIFF = rt.TH1D("phi_diff","phi_diff",40,0,1)
 
	for file_path in file_paths:
		try: 
			print "working on file", file_path
        		file = rt.TFile.Open(file_path)
        		# open tree file
        		tree = file.Get("demo/tree")
        		N = tree.GetEntries()
			print "nr of events:",N
        		#pt_diff, eta_diff, phi_diff = [],[],[]
        		
        		for i in xrange(N):
        		        if i % 100 == 0: print "Working on event " ,i
        		        if EarlyBreak > 0 and i>=EarlyBreak: break
        		        tree.GetEntry(i)
        		        for j in range(0,tree.nJets):
        		                jVector = rt.TLorentzVector()
        		                jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])
        		                quark_pt_eta_phi = Ana.SearchFor("quark", tree, jVector, dR, MomentumThreshold, BG)
        		                hadron_pt_eta_phi = Ana.SearchFor("hadron", tree, jVector, dR, MomentumThreshold, BG)
        		                if (quark_pt_eta_phi != False and hadron_pt_eta_phi != False):
        		                        pt_diff_temp = abs(quark_pt_eta_phi[0] - hadron_pt_eta_phi[0])/abs(quark_pt_eta_phi[0])
        		                        eta_diff_temp = abs(quark_pt_eta_phi[1] - hadron_pt_eta_phi[1])/abs(quark_pt_eta_phi[1])
        		                        phi_diff_temp = abs(quark_pt_eta_phi[2] - hadron_pt_eta_phi[2])/abs(quark_pt_eta_phi[2])
        		                        PT_DIFF.Fill(pt_diff_temp)
        		                        ETA_DIFF.Fill(eta_diff_temp)
        		                        PHI_DIFF.Fill(phi_diff_temp)
		except:
			print "skipped file because of error"
			continue

	tfile = rt.TFile("Thesis_Plots/root_files/QuarkvsHadron_{}.root".format(title),"recreate")
	PT_DIFF.Write()
        ETA_DIFF.Write()
        PHI_DIFF.Write()

def exclusive_tagged_jets_hist(signal_title, data, discriminant, discriminant_cut, CSV_cut,ran, nbins, Difference=False, mode="pT_jet",y_max=0,Save=False):
        """creates three histograms when given a signal dataset: -all jets tagged by CSV but not by the discriminant (L4-L1 for Difference = True, L4/L1 for Difference = False inside given dR), -all jets tagged
by the discriminant but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (see FeatureDict)"""
        title = signal_title+"_tagged_jets_vs_"+mode+"_exclusive"
        CSV_and_not_Discriminant = rt.TH1D(signal_title+"_CSV_not_Discriminant",title,nbins,ran[0],ran[1])
        CSV_and_not_Discriminant.SetLineColor(3)
        Discriminant_and_not_CSV = rt.TH1D(signal_title+"_Discriminant_and_not_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_not_CSV.SetLineColor(2)
        Discriminant_and_CSV = rt.TH1D(signal_title+"_Discriminant_and_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_CSV.SetLineColor(4)
	if Difference:
		title = title+"_delta"
	else:
		title = title+"_ratio"
        for particle in data:
                CSV_tag, Disc_tag = False, False
                if particle[1] >= CSV_cut:
                        CSV_tag = True
                if Difference:
                        L = particle[8]-particle[5]
                else:
                        if particle[17] != 0:
                                L = particle[20]/float(particle[17])
                        else:
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
	if y_max >0: CSV_and_not_Discriminant.SetMaximum(y_max)
        CSV_and_not_Discriminant.Draw()
        Discriminant_and_not_CSV.Draw("SAME")
	Discriminant_and_CSV.Draw("SAME")
        legend.Draw()
        if Save:
                canvas.SaveAs("Thesis_Plots/"+title+".png")
                Tfile= rt.TFile("Thesis_Plots/root_files/"+title+".root","recreate")
                Discriminant_and_not_CSV.Write()
                CSV_and_not_Discriminant.Write()
                Discriminant_and_CSV.Write()

def binned_exclusive_tagged_jets_hist(signal_title, data, discriminant, discriminant_cuts, CSV_cuts, bins, ran, nbins, Difference=False, mode="pT_jet", y_max=0,Save=False, Stacked=False, AllJets=None):
        """creates three histograms when given a signal dataset: -all jets tagged by CSV but not by the discriminant (L4-L1 for Difference = True, L4/L1 for Difference = False inside given dR), -all jets tagged
by the discriminant but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (see FeatureDict)"""
        title = signal_title+"_tagged_jets_vs_"+mode+"_exclusive"
        CSV_and_not_Discriminant = rt.TH1D(signal_title+"_CSV_not_Discriminant",title,nbins,ran[0],ran[1])
        CSV_and_not_Discriminant.SetLineColor(3)
        Discriminant_and_not_CSV = rt.TH1D(signal_title+"_Discriminant_and_not_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_not_CSV.SetLineColor(2)
        Discriminant_and_CSV = rt.TH1D(signal_title+"_Discriminant_and_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_CSV.SetLineColor(4)
	if Difference:
		title = title+"_delta"
	else:
		title = title+"_ratio"
        for particle in data:
		bin_number = FCM.bin_selection(particle,bins)
                if bin_number == -100: continue

                CSV_tag, Disc_tag = False, False
                if particle[1] >= CSV_cuts[bin_number]:
                        CSV_tag = True
                if Difference:
                        L = particle[8]-particle[5]
                else:
                        if particle[17] != 0:
                                L = particle[20]/float(particle[17])
                        else:
                               	L= -100
                if L >= discriminant_cuts[bin_number]:
                                Disc_tag = True
                if Disc_tag and not CSV_tag: Discriminant_and_not_CSV.Fill(particle[FeatureDict[mode]])
                if CSV_tag and not Disc_tag: CSV_and_not_Discriminant.Fill(particle[FeatureDict[mode]])
                if CSV_tag and Disc_tag: Discriminant_and_CSV.Fill(particle[FeatureDict[mode]])
	
	if Stacked:
		assert AllJets != None, "need to give a histogram containing full amount of jets into AllJets parameter"
		Stacked_tagged_jets(title, discriminant, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max)
	else:
        	canvas = rt.TCanvas("canvas","canvas",600,600)
        	canvas.SetTitle(title)
        	rt.gStyle.SetOptStat(0)
        	legend = rt.TLegend(0.9,0.9,0.65,0.75)
        	legend.AddEntry(Discriminant_and_not_CSV,discriminant+" and not CSV")
        	legend.AddEntry(CSV_and_not_Discriminant,"CSV and not "+discriminant)
        	legend.AddEntry(Discriminant_and_CSV, discriminant+" and CSV")
        	Discriminant_and_not_CSV.GetXaxis().SetTitle("jet p_{T} (GeV)")
        	Discriminant_and_not_CSV.GetYaxis().SetTitle('# jets')
        	Discriminant_and_not_CSV.GetYaxis().SetTitleOffset(1.5)
		if y_max >0: CSV_and_not_Discriminant.SetMaximum(y_max)
        	CSV_and_not_Discriminant.Draw()
        	Discriminant_and_not_CSV.Draw("SAME")
		Discriminant_and_CSV.Draw("SAME")
        	legend.Draw()
               	canvas.SaveAs("Thesis_Plots/"+title+".png")
	Tfile= rt.TFile("Thesis_Plots/root_files/"+title+".root","recreate")
	Discriminant_and_not_CSV.Write()
	CSV_and_not_Discriminant.Write()
	Discriminant_and_CSV.Write()

def ANN_predict(particle, model, model_type):
	"""takes the data of one particle (premade by Finalclustermatcher(), a keras model and a string describing its type and returns the b-tag value of that particle estimated by the given model"""
	hit_data = np.delete(particle[0:25],[0,1,2,3,4])
	conv_data = np.reshape(hit_data,(-1,5,4,1))

	if particle[17] != 0:
		L4_L1 = particle[20]/particle[17]
		L2_L1 = particle[18]/particle[17]
	else:
		L4_L1 = 10
		L2_L1 = 10
	if particle[18] != 0:
		L3_L2 = particle[19]/particle[18]
	else:
		L3_L2 = 10
	if particle[19] != 0:
		L4_L3 = particle[20]/particle[19]
	else:	
		L4_L3 = 10
		
	disc_data = np.array([[particle[8]-particle[5],L4_L1,L2_L1,L3_L2,L4_L3]])

	if model_type == "simple":
		input_data = hit_data
	elif model_type == "convolution":
		input_data = conv_data
	elif model_type == "functional":
		input_data = [conv_data,disc_data]
	elif model_type == "functional_pT":
		input_Data = [conv_data,disc_data,particle[3]]
	elif model_type == "functional_pV":
		input_Data = [conv_data,disc_data,particle[25]]

	return model.predict(input_data)[0]

def ANN_discriminants(x_data):
	"""helper function that takes the 5-cones hits data and returns an array containing different discriminants for the functional ANN"""
        L1_d = x_data[:,0]
        L4_d = x_data[:,3]
        L1_r = x_data[:,12]
        L2_r = x_data[:,13]
        L3_r = x_data[:,14]
        L4_r = x_data[:,15]
        L2_L1 = L2_r/L1_r
        L3_L2 = L3_r/L2_r
        L4_L3 = L4_r/L3_r
        L4_L1 = L4_r/L1_r
        L2_L1[np.isnan(L2_L1)]=1
        L3_L2[np.isnan(L3_L2)]=1
        L4_L3[np.isnan(L4_L3)]=1
        L4_L1[np.isnan(L4_L1)]=1
        L2_L1[L2_L1>100]=10
        L3_L2[L3_L2>100]=10
        L4_L3[L4_L3>100]=10
        L4_L1[L4_L1>100]=10
        return np.vstack((L4_d-L1_d,L4_L1,L2_L1,L3_L2,L4_L3)).transpose()

def ANN_functional_shape(x_data):
	"""takes the 5-cones hits data and puts it into the right input shape for the functional ANN model"""
	x_data_Li_Lj= ANN_discriminants(x_data)
	x_data_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
	return [x_data_Li, x_data_Li_Lj]

def ANN_clusterdata_to_x_pT_PV_CSV_label(signal_data,bg_data,shuffle=False):
	"""takes the data preprocessed by Finalclustermatcher() and returns a tuple of separated data necessary for any ANN model"""
	Signal_data = signal_data.copy()
	Signal_data[:,0]=1
	BG_data = bg_data.copy()
	BG_data[:,0] = 0	
	Data_sample = np.delete(np.vstack((Signal_data,BG_data)),[2,4],axis=1)
	if shuffle: np.random.shuffle(Data_sample)
	x_data = Data_sample[:,3:]
	y_data = Data_sample[:,0]
	CSV = Data_sample[:,1]
	pT = Data_sample[:,2]
	PV = Data_sample[:,-1]
	return (x_data, pT, PV, CSV, y_data)

def ANN_x_pT_CSV_label_to_clusterdata(x_data, pT, CSV, y_data):
	"""inverse of the above function: takes data necessary to the ANN and reverts it to the same shape as the data made by the Finalclustermatcher(). nEvent, hadron-pT and decayvx get lost"""
	blank = np.zeros(shape=(x_data.shape[0],1))
	Data_sample = np.hstack((y_data.reshape(-1,1),CSV.reshape(-1,1),blank,pT.reshape(-1,1),blank,x_data))
	signal = Data_sample[Data_sample[:,0]==1]
	bg = Data_sample[Data_sample[:,0]==0]
	return (signal, bg)

def ANN_bin_selection(pT,bins):
	"""numpy array version of the bin_selection() function: takes an array of all pT values and the pT-bins and returns an array of same length as pT with the corresponing bins index at each entry. pT values outside the bins are labeled with -100"""
	bin_numbers = np.zeros(len(pT))
	for n in range(len(bins)-1):
		bin_numbers += (n+100)*(pT>bins[n])*(pT<bins[n+1])
	bin_numbers -=100
	return bin_numbers.astype(int)

def ANN_Make_Binned_ROC_histograms(title,model, x_data, pT, CSV, bins, PU_range='full',addFeature=False):
	"""makes binned ROC histograms for an ANN. takes as input a keras model, the necessary ANN data, pT, CSV and the desired pT-binning"""
        nbins = 60

        ANN_hist_list = []
        CSV_hist_list = []
        for bin_ in range(len(bins)-1):
                ANN_hist_list.append(rt.TH1D("ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"ANN_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))
                CSV_hist_list.append(rt.TH1D("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),"CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1]),nbins,0,1))

	if addFeature == False:
		pred_y = model.predict(ANN_functional_shape(x_data))
	elif addFeature == "pT":
		pred_y = model.predict(ANN_functional_shape(x_data)+[pT/200])
	elif addFeature == "PV":
		assert x_data.shape[1] == 21, "wrong x_data shape: PV cannot be found"
		pred_y = model.predict(ANN_functional_shape(x_data)+[x_data[:,-1]/10.])
	else:
		print "invalid feature selection"
		return None
	bin_numbers = ANN_bin_selection(pT,bins)

        for n,particle in enumerate(x_data):
                if PU_range != 'full':
                        if particle[-1]<PU_range[0] or particle[-1]>PU_range[1]: continue
                if bin_numbers[n] == -100: continue
                ANN_hist_list[int(bin_numbers[n])].Fill(pred_y[n])
                CSV_hist_list[int(bin_numbers[n])].Fill(CSV[n])

        tfile = rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(title),"recreate")
        for hist in ANN_hist_list:
                hist.Write()
        for hist in CSV_hist_list:
                hist.Write()
        print "saved histograms in Thesis_Plots/root_files/{}_histograms.root".format(title)

def ANN_binned_tagged_jets_hist(datalist, model, discriminant_cuts, CSV_cuts, bins, nbins, mode="pT_jet",Save=False,addFeature=False):
        """creates a histogram for each dataset given as list of tuples (x_data, pT, CSV, title, range) of all the jets that were b-tagged by passing a given a list of cut values corresponding to the given pT-bins for CSV and ANN versus a feature (currently only working for jet-pT) given as string to 'mode'. The histograms are saved to a root file for further use."""
        title = "binned_tagged_jets_vs_"+mode
	discriminant = "ANN"
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
		datatitle = data[3]
                print "working on",datatitle
                ran = data[4]
		CSV = data[2]
		pT = data[1]
		x_data = data[0]
                AllJetsHistlist.append(rt.TH1D(datatitle+"_AllJets",datatitle+"_"+title,nbins,ran[0],ran[1]))
                AllJetsHistlist[n].SetLineColor(4)
                CSVHistlist.append(rt.TH1D(datatitle+"_CSV",datatitle+"_"+title,nbins,ran[0],ran[1]))
                CSVHistlist[n].SetLineColor(3)
                DiscriminantHistlist.append(rt.TH1D(datatitle+"_Discriminant",datatitle+"_"+title,nbins,ran[0],ran[1]))
                DiscriminantHistlist[n].SetLineColor(2)
	
		if addFeature == False:
			pred_y = model.predict(ANN_functional_shape(x_data))
		elif addFeature == "pT":
			pred_y = model.predict(ANN_functional_shape(x_data)+[pT/200])
		elif addFeature == "PV":
			assert x_data.shape[1] == 21, "wrong x_data format: PV cannot be found"
			pred_y = model.predict(ANN_functional_shape(x_data)+[x_data[:,-1]/10.])
		else:
			print "invalid feature input"
			return None
		bin_numbers = ANN_bin_selection(pT,bins)

	        for i,pT_value in enumerate(pT):
	                if bin_numbers[i] == -100: continue
			AllJetsHistlist[n].Fill(pT_value)
	                if pred_y[i] >= discriminant_cuts[bin_numbers[i]]: DiscriminantHistlist[n].Fill(pT_value)
	                if CSV[i] >= CSV_cuts[bin_numbers[i]]: CSVHistlist[n].Fill(pT_value)

        canvaslist = []
        legendlist = []
        Tfilelist = []
        for n,data in enumerate(datalist):
		datatitle = data[3]
                canvaslist.append(rt.TCanvas(datatitle+"_canvas","canvas",600,600))
                canvaslist[n].SetTitle(datatitle+"_"+title)
                rt.gStyle.SetOptStat(0)
                legendlist.append(rt.TLegend(0.9,0.9,0.65,0.75))
                legendlist[n].AddEntry(AllJetsHistlist[n], "All jets")
                legendlist[n].AddEntry(CSVHistlist[n], "CSV")
                legendlist[n].AddEntry(DiscriminantHistlist[n], discriminant)
                AllJetsHistlist[n].GetXaxis().SetTitle(mode)
                AllJetsHistlist[n].GetYaxis().SetTitle('# jets')
                AllJetsHistlist[n].GetYaxis().SetTitleOffset(1.5)
                #AllJetsHistlist[n].Draw()
                #CSVHistlist[n].Draw("SAME")
                #DiscriminantHistlist[n].Draw("SAME")
                #legendlist[n].Draw()
                if Save:
                        #canvaslist[n].SaveAs(title+"_"+datatitle+discriminant+".png")
                        Tfilelist.append(rt.TFile("Thesis_Plots/root_files/"+title+"_"+datatitle+discriminant+".root","recreate"))
                        print "saved histogram as Thesis_Plots/root_files/"+title+"_"+datatitle+discriminant+".root"
                        AllJetsHistlist[n].Write()
                        CSVHistlist[n].Write()
                        DiscriminantHistlist[n].Write()

def ANN_exclusive_tagged_jets_hist(signal_title, model, x_data, pT, CSV, discriminant_cuts, CSV_cuts, bins, ran, nbins, mode="pT_jet",y_max=0,Save=False, addFeature=False, DrawTitle=False, Stacked=False, AllJets=None):
        """creates three histograms when given a signal dataset: -all jets tagged by CSV but not by ANN, -all jets tagged
by ANN but not by CSV, -all jets tagged by both. the amount of jets is plotted versus the feature given as string to 'mode' (currently only jet-pT working) (see FeatureDict)"""
        title = signal_title+"_tagged_jets_vs_"+mode+"_exclusive"
	discriminant = "ANN"
        CSV_and_not_Discriminant = rt.TH1D(signal_title+"_CSV_not_Discriminant",title,nbins,ran[0],ran[1])
        CSV_and_not_Discriminant.SetLineColor(3)
        Discriminant_and_not_CSV = rt.TH1D(signal_title+"_Discriminant_and_not_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_not_CSV.SetLineColor(2)
        Discriminant_and_CSV = rt.TH1D(signal_title+"_Discriminant_and_CSV",title,nbins,ran[0],ran[1])
        Discriminant_and_CSV.SetLineColor(4)
       
	if addFeature == False:
		pred_y = model.predict(ANN_functional_shape(x_data))
	elif addFeature == "pT":
		pred_y = model.predict(ANN_functional_shape(x_data)+[pT/200])
	elif addFeature == "PV":
		assert x_data.shape[1] == 21, "wrong x_data format: PV cannot be found"
		pred_y = model.predict(ANN_functional_shape(x_data)+[x_data[:,-1]/10])
	else:
		print "invalid feature input"
		return None

	bin_numbers = ANN_bin_selection(pT,bins)

	for i,pT_value in enumerate(pT):
	        if bin_numbers[i] == -100: continue
                CSV_tag, Disc_tag = False, False
                if CSV[i] >= CSV_cuts[bin_numbers[i]]:
                        CSV_tag = True
                
		if pred_y[i] >= discriminant_cuts[bin_numbers[i]]:
                                Disc_tag = True
                if Disc_tag and not CSV_tag: Discriminant_and_not_CSV.Fill(pT_value)
                if CSV_tag and not Disc_tag: CSV_and_not_Discriminant.Fill(pT_value)
                if CSV_tag and Disc_tag: Discriminant_and_CSV.Fill(pT_value)

	if Stacked:
		assert AllJets != None, "need to give a histogram containing full amount of jets into AllJets parameter"
		Stacked_tagged_jets(title, discriminant, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max)
	else:
        	canvas = rt.TCanvas("canvas","canvas",600,600)
        	if DrawTitle: 
			canvas.SetTitle(title)
		else:
			rt.gStyle.SetOptTitle(0)
        	rt.gStyle.SetOptStat(0)
        	legend = rt.TLegend(0.9,0.9,0.65,0.75)
        	legend.AddEntry(Discriminant_and_not_CSV,discriminant+" and not CSV")
        	legend.AddEntry(CSV_and_not_Discriminant,"CSV and not "+discriminant)
        	legend.AddEntry(Discriminant_and_CSV, discriminant+" and CSV")
        	CSV_and_not_Discriminant.GetXaxis().SetTitle("jet p_{T} (GeV)")
        	CSV_and_not_Discriminant.GetYaxis().SetTitle('# jets')
        	CSV_and_not_Discriminant.GetYaxis().SetTitleOffset(1.5)
		if y_max >0: CSV_and_not_Discriminant.SetMaximum(y_max)
        	CSV_and_not_Discriminant.Draw()
        	Discriminant_and_not_CSV.Draw("SAME")
		Discriminant_and_CSV.Draw("SAME")
        	legend.Draw()
                canvas.SaveAs("Thesis_Plots/"+title+".png")
        Tfile= rt.TFile("Thesis_Plots/root_files/"+title+".root","recreate")
        Discriminant_and_not_CSV.Write()
        CSV_and_not_Discriminant.Write()
        Discriminant_and_CSV.Write()
	print "saved as Thesis_Plots/root_files/"+title+".root"

def Make_Binned_ANN_ROC_Curves(title,Signal_title,Background_title,bins,log=False):
	"""accesses the files made by ANN_Make_Binned_ROC_histograms() directly to produce a ROC curve for ANN and CSV in the desired pT-bins. A log representation can be turned on"""
        #hsv = plt.get_cmap('hsv')
        #color = hsv(np.linspace(0,1.0,len(bins)-1))
        #color = ['b', 'g', 'r', 'c', 'm', 'y']
        if len(bins)<=6:
                color = ['red','green','blue','orange','brown']
        else:
                color = ['deepskyblue','rosybrown','olivedrab','royalblue','firebrick','chartreuse','navy','red','darkorchid','lightseagreen','mediumvioletred','blue']
        nbins = 60
	dis_string = "ANN_"

        Signal_file = rt.TFile("Thesis_Plots/root_files/{}_ANN_histograms.root".format(Signal_title),"READ")
        Background_file =   rt.TFile("Thesis_Plots/root_files/{}_ANN_histograms.root".format(Background_title),"READ")

        plt.figure("ROC")
        plt.clf()

        for bin_ in range(len(bins)-1):
                Dis_Signal_Eff = FCM.Get_ROC_Efficiencies(Signal_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,0)
                Dis_BG_Eff = FCM.Get_ROC_Efficiencies(Background_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),nbins,0)
                CSV_Signal_Eff = FCM.Get_ROC_Efficiencies(Signal_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),ratio_bins,0)
                CSV_BG_Eff = FCM.Get_ROC_Efficiencies(Background_file.Get("CSV_"+str(bins[bin_])+"_"+str(bins[bin_+1])),(0,1),ratio_bins,0)
                if log:
                        plt.semilogy(Dis_Signal_Eff,Dis_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"_"+str(bins[bin_+1]))
                        plt.semilogy(CSV_Signal_Eff,CSV_BG_Eff, color = color[bin_],linestyle = '--',)

                else:
                        plt.plot(Dis_Signal_Eff,1-Dis_BG_Eff, color = color[bin_], linestyle = '-',label=str(bins[bin_])+"_"+str(bins[bin_+1]))
                        plt.plot(CSV_Signal_Eff,1-CSV_BG_Eff, color = color[bin_],linestyle = '--',)

        if log:
		if diff:
			plt.semilogy([0,0],[0,0],'k-',label = 'L4-L1')
		else:
                	plt.semilogy([0,0],[0,0],'k-',label = 'L4/L1')
                plt.semilogy([0,0],[0,0],'k-.',label = 'CSV')
                plt.semilogy([0,1],[0.1,0.1],'k:')
                plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"$\epsilon$_background")
                plt.legend(loc=4)
        else:
		if diff:
			plt.plot([0,0],[0,0],'k-',label = 'L4-L1')
		else:
                	plt.plot([0,0],[0,0],'k-',label = 'L4/L1')
                plt.plot([0,0],[0,0],'k-.',label = 'CSV')
                #plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
                plt.plot([0,1],[0.9,0.9],'k:')
                plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"1-$\epsilon$_background")
                plt.legend(loc=3)
        #plt.title(title+"_ROC-Curves")

        plt.savefig("Thesis_Plots/{}_ROC_Curves.png".format(title))
        print "saved as Thesis_Plots/{}_ROC_Curves.png".format(title)

def ANN_efficiency_vs_PU(title, x_data, pT, CSV, model, ANN_Cuts, Ratio_Cuts, CSV_Cuts, bins, y_max, pT_Cut=200, BG=False, DrawTitle=False):
	"""draws the efficiency vs number of PV for ANN, L4/L1 and CSV when given the ANN data and the cuts corresponding to each discriminant"""
        assert x_data.shape[1]==21, "x_data does not contain PV. Make sure it is made from a PU sample and has shape (x, 21)."
	assert x_data.shape[0] == len(pT) == len(CSV), "data inputs need to have the same length"
	assert len(ANN_Cuts) == len(Ratio_Cuts) == len(CSV_Cuts) == len(bins)-1, "cuts need to have the same length and be compatible with amount of bins"

        ran = (0,80)
        nbins = 80
        import array
	if BG:
		bins_ = array.array('d',[0.0, 11.0]+range(19,41,8)+[42.0,  52.0, 80])
	else:
        	bins_ = array.array('d',[0.0, 11.0]+range(15,41,4)+[42.0, 52.0, 58.0, 65.0, 80])

	if pT_Cut >= 1200:
		bins_ = array.array('d',[0.0, 20.0, 40.0, 80.0])


        #make histograms of efficiency vs PU
        AllJets_Hist = rt.TH1D("AllJets","AllJets",nbins,ran[0],ran[1])
        ANN_Hist = rt.TH1D("ANN","ANN",nbins,ran[0],ran[1])
        Ratio_Hist = rt.TH1D("Ratio","Ratio",nbins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",nbins,ran[0],ran[1])

	AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        ANN_Hist = ANN_Hist.Rebin(len(bins_)-1,"ANN",bins_)
        Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)
 
	pred_y = model.predict(ANN_functional_shape(x_data))
	bin_numbers = ANN_bin_selection(pT,bins)
	
	for i,pT_value in enumerate(pT):
			if pT_value < pT_Cut: continue
	                if bin_numbers[i] == -100: continue
			AllJets_Hist.Fill(x_data[i,-1])
			if CSV[i] >= CSV_Cuts[bin_numbers[i]]: CSV_Hist.Fill(x_data[i,-1])
	                if pred_y[i] >= ANN_Cuts[bin_numbers[i]]: ANN_Hist.Fill(x_data[i,-1])
			if x_data[i,12] != 0:
				L_R = x_data[i,15]/float(x_data[i,12])
				if L_R >= Ratio_Cuts[bin_numbers[i]]: Ratio_Hist.Fill(x_data[i,-1])
	                
	'''		
        AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        ANN_Hist = ANN_Hist.Rebin(len(bins_)-1,"ANN",bins_)
        Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)
	'''
        #Make Graphs and draw them
        canvas = rt.TCanvas('canvas','canvas',600,600)
	if DrawTitle == False: rt.gStyle.SetOptTitle(0)
        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        ANN_Graph = rt.TGraphAsymmErrors()
        Ratio_Graph = rt.TGraphAsymmErrors()
        CSV_Graph = rt.TGraphAsymmErrors()
        if DrawTitle: Ratio_Graph.SetTitle(title+"_vs_PU_pT{}{}".format('jet',pT_Cut))
        ANN_Graph.Divide(ANN_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        ANN_Graph.SetLineColor(3)
        Ratio_Graph.SetLineColor(2)
        CSV_Graph.SetLineColor(4)
        legend.AddEntry(ANN_Graph, "ANN", "LEP")
        legend.AddEntry(Ratio_Graph, "L4/L1", "LEP")
        legend.AddEntry(CSV_Graph, "CSV", "LEP")
        Ratio_Graph.GetXaxis().SetTitle("#PV")
        if BG:
		Ratio_Graph.GetYaxis().SetTitle('mistag rate')
	else:	
		Ratio_Graph.GetYaxis().SetTitle('efficiency')
        Ratio_Graph.GetYaxis().SetTitleOffset(1.5)
	Ratio_Graph.SetMinimum(0.)
        Ratio_Graph.SetMaximum(y_max)
        Ratio_Graph.Draw()
        ANN_Graph.Draw("SAME")
        CSV_Graph.Draw("SAME")
        legend.Draw()
        canvas.SaveAs('Thesis_Plots/'+title+"_vs_PU_pT{}{}.png".format('jet',pT_Cut))


def ANN_efficiency_vs_PU_pT_PV(title, x_data, pT, CSV, model_noPT, model_withPT, model_withPV, ANN_noPT_Cuts, ANN_withPT_Cuts, ANN_withPV_Cuts, Ratio_Cuts, CSV_Cuts, bins, y_max, pT_Cut=200, BG=False, DrawTitle=False, LargeLegend=False):
	"""draws the efficiency vs number of PV for ANN, L4/L1 and CSV when given the ANN data and the cuts corresponding to each discriminant"""
        assert x_data.shape[1]==21, "x_data does not contain PV. Make sure it is made from a PU sample and has shape (x, 21)."
	assert x_data.shape[0] == len(pT) == len(CSV), "data inputs need to have the same length"
	assert len(ANN_noPT_Cuts) == len(ANN_withPT_Cuts) == len(ANN_withPV_Cuts) == len(Ratio_Cuts) == len(CSV_Cuts) == len(bins)-1, "cuts need to have the same length and be compatible with amount of bins"

        ran = (0,80)
        nbins = 80
        import array
	if BG:
		bins_ = array.array('d',[0.0, 11.0]+range(19,41,8)+[42.0,  52.0, 80])
	else:
        	bins_ = array.array('d',[0.0, 11.0]+range(15,41,4)+[42.0, 52.0, 58.0, 65.0, 80])

	if pT_Cut >= 1200:
		bins_ = array.array('d',[0.0, 20.0, 40.0, 80.0])


        #make histograms of efficiency vs PU
        AllJets_Hist = rt.TH1D("AllJets","AllJets",nbins,ran[0],ran[1])
        ANN_noPT_Hist = rt.TH1D("ANN_noPT","ANN_noPT",nbins,ran[0],ran[1])
	ANN_withPT_Hist = rt.TH1D("ANN_withPT","ANN_withPT",nbins,ran[0],ran[1])
	ANN_withPV_Hist = rt.TH1D("ANN_withPV","ANN_withPV",nbins,ran[0],ran[1])
        Ratio_Hist = rt.TH1D("Ratio","Ratio",nbins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",nbins,ran[0],ran[1])

	AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        ANN_noPT_Hist = ANN_noPT_Hist.Rebin(len(bins_)-1,"ANN_noPT",bins_)
	ANN_withPT_Hist = ANN_withPT_Hist.Rebin(len(bins_)-1,"ANN_withPT",bins_)
	ANN_withPV_Hist = ANN_withPV_Hist.Rebin(len(bins_)-1,"ANN_withPV",bins_)
        Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)
 
	pred_y_noPT = model_noPT.predict(ANN_functional_shape(x_data))
	pred_y_withPT = model_withPT.predict(ANN_functional_shape(x_data)+[pT/200.])
	pred_y_withPV = model_withPV.predict(ANN_functional_shape(x_data)+[x_data[:,-1]/10.])

	bin_numbers = ANN_bin_selection(pT,bins)
	
	for i,pT_value in enumerate(pT):
			if pT_value < pT_Cut: continue
	                if bin_numbers[i] == -100: continue
			AllJets_Hist.Fill(x_data[i,-1])
			if CSV[i] >= CSV_Cuts[bin_numbers[i]]: CSV_Hist.Fill(x_data[i,-1])
	                if pred_y_noPT[i] >= ANN_noPT_Cuts[bin_numbers[i]]: ANN_noPT_Hist.Fill(x_data[i,-1])
			if pred_y_withPT[i] >= ANN_withPT_Cuts[bin_numbers[i]]: ANN_withPT_Hist.Fill(x_data[i,-1])
			if pred_y_withPV[i] >= ANN_withPV_Cuts[bin_numbers[i]]: ANN_withPV_Hist.Fill(x_data[i,-1])

			if x_data[i,12] != 0:
				L_R = x_data[i,15]/float(x_data[i,12])
				if L_R >= Ratio_Cuts[bin_numbers[i]]: Ratio_Hist.Fill(x_data[i,-1])
        
	#Make Graphs and draw them
        canvas = rt.TCanvas('canvas','canvas',600,600)
	if DrawTitle == False: rt.gStyle.SetOptTitle(0)
	if LargeLegend:
		legend = rt.TLegend(0.1,0.9,0.4,0.7)
	else:
        	legend = rt.TLegend(0.1,0.9,0.35,0.75)
        ANN_noPT_Graph = rt.TGraphAsymmErrors()
	ANN_withPT_Graph = rt.TGraphAsymmErrors()
	ANN_withPV_Graph = rt.TGraphAsymmErrors()
        Ratio_Graph = rt.TGraphAsymmErrors()
        CSV_Graph = rt.TGraphAsymmErrors()
        if DrawTitle: Ratio_Graph.SetTitle(title+"_vs_PU_pT{}{}".format('jet',pT_Cut))
        ANN_noPT_Graph.Divide(ANN_noPT_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	ANN_withPT_Graph.Divide(ANN_withPT_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
	ANN_withPV_Graph.Divide(ANN_withPV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        ANN_noPT_Graph.SetLineColor(3)
        ANN_withPT_Graph.SetLineColor(6)
	ANN_withPV_Graph.SetLineColor(7)
	Ratio_Graph.SetLineColor(2)
        CSV_Graph.SetLineColor(4)
        #legend.AddEntry(ANN_noPT_Graph, "ANN without p_{T}/PV", "LEP")
	legend.AddEntry(ANN_noPT_Graph, "ANN without p_{T}", "LEP")
        legend.AddEntry(ANN_withPT_Graph, "ANN with p_{T}", "LEP")
	#legend.AddEntry(ANN_withPV_Graph, "ANN with PV", "LEP")
	legend.AddEntry(Ratio_Graph, "L4/L1", "LEP")
        legend.AddEntry(CSV_Graph, "CSV", "LEP")
        Ratio_Graph.GetXaxis().SetTitle("#PV")
        if BG:
		Ratio_Graph.GetYaxis().SetTitle('mistag rate')
	else:	
		Ratio_Graph.GetYaxis().SetTitle('efficiency')
        Ratio_Graph.GetYaxis().SetTitleOffset(1.5)
	Ratio_Graph.SetMinimum(0.)
        Ratio_Graph.SetMaximum(y_max)
        Ratio_Graph.Draw()
        ANN_noPT_Graph.Draw("SAME")
	ANN_withPT_Graph.Draw("SAME")
	#ANN_withPV_Graph.Draw("SAME")
        CSV_Graph.Draw("SAME")
        legend.Draw()
        canvas.SaveAs('Thesis_Plots/'+title+"_vs_PU_pT{}{}.png".format('jet',pT_Cut))

def Stacked_tagged_jets(title, Discriminant_Name, AllJets, Discriminant_and_not_CSV, CSV_and_not_Discriminant, Discriminant_and_CSV, y_max):
	
	AllJets_hist = AllJets.Clone()
	Discriminant_and_not_CSV_hist = Discriminant_and_not_CSV.Clone()
	CSV_and_not_Discriminant_hist = CSV_and_not_Discriminant.Clone()
	Discriminant_and_CSV_hist = Discriminant_and_CSV.Clone()

	stack = rt.THStack("stack", "stack")
	AllJets_hist.SetLineColor(1)
	Discriminant_and_not_CSV_hist.SetFillColor(2)
	CSV_and_not_Discriminant_hist.SetFillColor(3)
	Discriminant_and_CSV_hist.SetFillColor(4)
	stack.Add(Discriminant_and_CSV_hist)
	stack.Add(CSV_and_not_Discriminant_hist)
	stack.Add(Discriminant_and_not_CSV_hist)

	canvas = rt.TCanvas("canvas","canvas",600,600)
	rt.gStyle.SetOptTitle(0)
        rt.gStyle.SetOptStat(0)
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
	legend.AddEntry(AllJets_hist,"all jets")
        legend.AddEntry(Discriminant_and_not_CSV_hist,Discriminant_Name+" and not CSV")
        legend.AddEntry(CSV_and_not_Discriminant_hist,"CSV and not "+Discriminant_Name)
        legend.AddEntry(Discriminant_and_CSV_hist, Discriminant_Name+" and CSV")
        AllJets_hist.GetXaxis().SetTitle("jet p_{T} (GeV)")
        AllJets_hist.GetYaxis().SetTitle('# jets')
        AllJets_hist.GetYaxis().SetTitleOffset(1.5)
	if y_max >0: AllJets_hist.SetMaximum(y_max)
        AllJets_hist.Draw()
	stack.Draw("SAME")
	rt.gPad.RedrawAxis()
        legend.Draw()
        canvas.SaveAs("Thesis_Plots/stacked_"+title+".png")

def Relative_Gains(Discriminant_not_CSV, CSV_not_Discriminant, Discriminant_and_CSV):
	ran = (0,2500)
	res = 60
	gain_list = []
	errors = []
	CSV = CSV_not_Discriminant.Clone()
	CSV.Add(Discriminant_and_CSV)
	thresholds = np.linspace(ran[0],ran[1],60)
	for threshold in thresholds:
		N_Disc_excl = Discriminant_not_CSV.Integral(Discriminant_not_CSV.GetXaxis().FindBin(threshold),Discriminant_not_CSV.GetXaxis().FindBin(ran[1]))
		N_CSV = CSV.Integral(CSV.GetXaxis().FindBin(threshold),CSV.GetXaxis().FindBin(ran[1]))
		if N_CSV != 0:
			gain_list.append(N_Disc_excl/N_CSV)
			errors.append((np.sqrt(N_Disc_excl)/N_CSV) * np.sqrt(1+ N_Disc_excl/N_CSV))
		else:
			gain_list.append(np.nan)
			errors.append(np.nan)
	return (np.array(gain_list), np.array(errors), thresholds)

def Relative_Gain_Plots_binned():
	
	ratio_file = rt.TFile.Open("Thesis_Plots/root_files/Signal_binned_tagged_jets_vs_pT_jet_exclusive_ratio.root")
	delta_file = rt.TFile.Open("Thesis_Plots/root_files/Signal_binned_tagged_jets_vs_pT_jet_exclusive_delta.root")

	ratio_not_CSV = ratio_file.Get("Signal_binned_Discriminant_and_not_CSV")
	delta_not_CSV = delta_file.Get("Signal_binned_Discriminant_and_not_CSV")
	CSV_not_ratio = ratio_file.Get("Signal_binned_CSV_not_Discriminant")
	CSV_not_delta = delta_file.Get("Signal_binned_CSV_not_Discriminant")
	ratio_and_CSV = ratio_file.Get("Signal_binned_Discriminant_and_CSV")
	delta_and_CSV = delta_file.Get("Signal_binned_Discriminant_and_CSV")
	
	ratio_gains, ratio_errors, thresholds = Relative_Gains(ratio_not_CSV, CSV_not_ratio, ratio_and_CSV)
	delta_gains, delta_errors, thresholds = Relative_Gains(delta_not_CSV, CSV_not_delta, delta_and_CSV)
	
	plt.figure()
	plt.errorbar(thresholds, ratio_gains, yerr=ratio_errors, fmt='r', label=r"L4/L1")
	plt.errorbar(thresholds, delta_gains, yerr=delta_errors, fmt='g', label=r"L4-L1")
	plt.xlim(200,2000)
	plt.ylim(0,1.2)
	plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
	plt.legend(loc=2)
	plt.savefig("Thesis_Plots/Relative_Gain_Plots_binned.png")
	print "saved figure as Thesis_Plots/Relative_Gain_Plots_binned.png"
	plt.show()

	
def Relative_Gain_Plots_ANN():
	
	ANN_noPT_file = rt.TFile.Open("Thesis_Plots/root_files/Signal_ANN_noPT_tagged_jets_vs_pT_jet_exclusive.root")
	ANN_withPT_file = rt.TFile.Open("Thesis_Plots/root_files/Signal_ANN_withPT_tagged_jets_vs_pT_jet_exclusive.root")
	L4_L1_file = rt.TFile.Open("Thesis_Plots/root_files/Signal_compare_tagged_jets_vs_pT_jet_exclusive_ratio.root")

	noPT_Discriminant_not_CSV = ANN_noPT_file.Get("Signal_ANN_noPT_Discriminant_and_not_CSV")
	noPT_CSV_not_Discriminant = ANN_noPT_file.Get("Signal_ANN_noPT_CSV_not_Discriminant")
	noPT_Discrimant_and_CSV = ANN_noPT_file.Get("Signal_ANN_noPT_Discriminant_and_CSV")
	withPT_Discriminant_not_CSV = ANN_withPT_file.Get("Signal_ANN_withPT_Discriminant_and_not_CSV")
	withPT_CSV_not_Discriminant = ANN_withPT_file.Get("Signal_ANN_withPT_CSV_not_Discriminant")
	withPT_Discrimant_and_CSV = ANN_withPT_file.Get("Signal_ANN_withPT_Discriminant_and_CSV")
	L4_L1_Discriminant_not_CSV = L4_L1_file.Get("Signal_compare_Discriminant_and_not_CSV")
	L4_L1_CSV_not_Discriminant = L4_L1_file.Get("Signal_compare_CSV_not_Discriminant")
	L4_L1_Discrimant_and_CSV = L4_L1_file.Get("Signal_compare_Discriminant_and_CSV")

	ANN_noPT_gains, ANN_noPT_errors, thresholds = Relative_Gains(noPT_Discriminant_not_CSV, noPT_CSV_not_Discriminant, noPT_Discrimant_and_CSV)
	ANN_withPT_gains, ANN_withPT_errors, thresholds = Relative_Gains(withPT_Discriminant_not_CSV, withPT_CSV_not_Discriminant, withPT_Discrimant_and_CSV)	
	L4_L1_gains, L4_L1_errors, thresholds = Relative_Gains(L4_L1_Discriminant_not_CSV, L4_L1_CSV_not_Discriminant, L4_L1_Discrimant_and_CSV)
	
	plt.figure()
	plt.errorbar(thresholds, ANN_noPT_gains, yerr=ANN_noPT_errors, fmt='g', label=r"ANN without $p_T$")
	plt.errorbar(thresholds, ANN_withPT_gains, yerr=ANN_withPT_errors, fmt='magenta', label=r"ANN with $p_T$")
	plt.errorbar(thresholds, L4_L1_gains, yerr=L4_L1_errors, fmt='r', label=r"L4/L1")
	plt.xlim(200,2000)
	plt.ylim(0,3.5)
	plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
	plt.legend(loc=2)
	plt.savefig("Thesis_Plots/Relative_Gain_Plots.png")
	print "saved figure as Thesis_Plots/Relative_Gain_Plots.png"
	plt.show()

def Draw_dR_dist_Histograms(title, L1_hist,L2_hist,L3_hist,L4_hist):
        canvas = rt.TCanvas('canvas','canvas',600,600)
	ran=(0,0.2)
	xlabel = '#DeltaR'
	ylabel = '(a.u.)'

	rt.gStyle.SetOptTitle(0)
	
        rt.gStyle.SetOptStat(0)#something is wrong with this
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
	L1 = L1_hist.Clone()
	L2 = L2_hist.Clone()	
	L3 = L3_hist.Clone()
	L4 = L4_hist.Clone()
	L1.SetLineColor(2)
	L2.SetLineColor(3)
	L3.SetLineColor(4)
	L4.SetLineColor(6)
	legend.AddEntry(L1,'Layer 1')
	legend.AddEntry(L2,'Layer 2')
	legend.AddEntry(L3,'Layer 3')
	legend.AddEntry(L4,'Layer 4')
	L4.GetXaxis().SetTitle(xlabel)
        L4.GetYaxis().SetTitle(ylabel)
        L4.GetYaxis().SetTitleOffset(1.5)
        L4.DrawNormalized()
	L1.DrawNormalized("SAME")
	L2.DrawNormalized("SAME")
	L3.DrawNormalized("SAME")
	legend.Draw()
        canvas.SaveAs("Thesis_Plots/dR_dist_"+title+".png")

def decayvx_plots(hist_2TeV, hist_4TeV, hist_BG):
	ran = (0,35)
	hist2TeV = hist_2TeV.Clone()
	hist4TeV = hist_4TeV.Clone()
	histBG = hist_BG.Clone()
	canvas = rt.TCanvas('canvas','canvas',600,600)
        rt.gStyle.SetOptTitle(0)
	legend = rt.TLegend(0.9,0.9,0.65,0.75)
	histBG.SetLineColor(2)
	hist2TeV.SetLineColor(3)
        hist4TeV.SetLineColor(4)
	legend.AddEntry(histBG, "background")
	legend.AddEntry(hist2TeV, "2TeV signal")
	legend.AddEntry(hist4TeV, "4TeV signal")
	hist2TeV.GetXaxis().SetTitle("decay vertex R (cm)")
	hist2TeV.GetYaxis().SetTitle("(a.u.)")
	hist2TeV.GetYaxis().SetTitleOffset(1.5)
	hist2TeV.DrawNormalized()
	hist4TeV.DrawNormalized("SAME")
	histBG.DrawNormalized("SAME")
	legend.Draw()
	canvas.SaveAs("Thesis_Plots/decay_vertex_R_Signal.png")


def DrawHistograms(Histograms, ran, title, xlabel, ylabel, Save=False,Normalize=True,DrawTitle=False, t_sleep=0):
        """Draws multiple histograms neatly on a canvas when given to the function as list of tuples (root histogram, title, colorindex(optional))."""
        canvas = rt.TCanvas('canvas','canvas',600,600)
	if DrawTitle: 
		canvas.SetTitle(title)
	else:
		rt.gStyle.SetOptTitle(0)
	
	histlist = []
        if len(Histograms) > 1:
                rt.gStyle.SetOptStat(0)#something is wrong with this
                legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for nr, Histogram in enumerate(Histograms):
		histlist.append(Histogram[0])
		if len(Histogram)>2:
			histlist[nr].SetLineColor(Histogram[2])
		else:
			if nr < 3:
               		 	histlist[nr].SetLineColor(nr+2)
			else:
				histlist[nr].SetLineColor(nr+3)
                if nr == 0:
			if DrawTitle: histlist[nr].SetTitle(title)
                        histlist[nr].GetXaxis().SetTitle(xlabel)
                        histlist[nr].GetYaxis().SetTitle(ylabel)
                        histlist[nr].GetYaxis().SetTitleOffset(1.5)
                        if Normalize:
                        	histlist[nr].DrawNormalized()
                        else:
                                histlist[nr].Draw()
                else:
                        if Normalize:
                                histlist[nr].DrawNormalized("SAME")
                        else:
                                histlist[nr].Draw("SAME")
                if len(Histograms)>1:
                        legend.AddEntry(histlist[nr],Histogram[1])
        if len(Histograms)>1: 
		#rt.gStyle.SetOptStat(0)#something is wrong with this
		legend.Draw()
        if Save: canvas.SaveAs("Thesis_Plots/"+title+".png")
        sleep(t_sleep)

if __name__ == '__main__':

	#load necessary data
	
	'''Data Strings'''
	
	#Signal_2TeV_noPU_String = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M2000_GENSIMDIGIRECO_v2.root"
	#Signal_4TeV_noPU_String = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M4000_GENSIMDIGIRECO_v2.root"
	Signal_2TeV_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_2000/btagHits_noPU/180628_091836/0000/flatTuple_{}.root' #18-33
        Signal_4TeV_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_noPU/180628_092048/0000/flatTuple_{}.root' #18-33
	Signal_2TeV_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_2000/btagHits_wPU/180604_113337/0000/flatTuple_{}.root' #1-17
	Signal_4TeV_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/ZprimeToBBbar_M_4000/btagHits_wPU/180604_071651/0000/flatTuple_{}.root' #1-16
	BG_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_wPVs/180502_130824/0000/flatTuple_{}.root' #1-47
	BG_PU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_wPU-v2/180528_101050/0000/flatTuple_{}.root' #1-500

	Signal_2TeV_noPU_String_list=[]
        Signal_4TeV_noPU_String_list=[]
        Signal_2TeV_PU_String_list=[]
	for i in range(18,34):
		Signal_2TeV_noPU_String_list.append(Signal_2TeV_noPU_String.format(i))
	for i in range(18,34):
		Signal_4TeV_noPU_String_list.append(Signal_4TeV_noPU_String.format(i))
	for i in range(1,18):
		Signal_2TeV_PU_String_list.append(Signal_2TeV_PU_String.format(i))
        Signal_4TeV_PU_String_list=[]
        for i in range(1,17):
		Signal_4TeV_PU_String_list.append(Signal_4TeV_PU_String.format(i))
	BG_noPU_String_list=[]
        for i in range(1,48):
		BG_noPU_String_list.append(BG_noPU_String.format(i))
	BG_PU_String_list=[]
	for i in range(1,501):
		BG_PU_String_list.append(BG_PU_String.format(i))
	
	
	#Signal_4TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_4TeV-Signal.npy')
        #Signal_2TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_2TeV-Signal.npy')
        #Signal_both_noPU = np.vstack((Signal_4TeV_noPU,Signal_2TeV_noPU))
        #Background_noPU = FCM.load_data('BG','matched_clusters/BG_noPU/',31, old_version=True)#31
        
	'''		
	Signal_4TeV_PU = FCM.load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
        Signal_2TeV_PU = FCM.load_data('2TeV-Signal_PU','matched_clusters/Signal_PU/',19)
        Signal_both_PU = np.vstack((Signal_4TeV_PU,Signal_2TeV_PU))
        Background_PU = FCM.load_data('BG_PU','matched_clusters/BG_PU/',499)	
	'''
	#Analysis of samples: mass?, pT, eta?, decayvx, CSV on quarks hadrons jets

	'''	
	#Efficient_Sample_Analysis("Signal_2TeV_jet_pT", Signal_2TeV_noPU, (0,2500), "jet_pT")
	#Efficient_Sample_Analysis("Signal_4TeV_jet_pT", Signal_4TeV_noPU, (0,2500), "jet_pT")
	#Efficient_Sample_Analysis("BG_jet_pT", Background_noPU, (0,2500), "jet_pT")
	tfile1 = rt.TFile.Open("Thesis_Plots/root_files/Signal_2TeV_jet_pT.root")
	jet_pT_hist_2TeV = tfile1.Get("jet_pT")
	tfile2 = rt.TFile.Open("Thesis_Plots/root_files/Signal_4TeV_jet_pT.root")
	jet_pT_hist_4TeV = tfile2.Get("jet_pT")
	tfile3 = rt.TFile.Open("Thesis_Plots/root_files/BG_jet_pT.root")
	jet_pT_hist_BG = tfile3.Get("jet_pT")
	Histograms = [(jet_pT_hist_BG,"background"),(jet_pT_hist_2TeV,"2TeV signal"),(jet_pT_hist_4TeV,"4TeV signal")]
	DrawHistograms(Histograms, (0,2500), "jet_pT", "jet p_{T} (GeV)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''
	'''
	#Efficient_Sample_Analysis("Signal_2TeV_decayvx", Signal_2TeV_noPU, (0,35), "decayvx")
	#Efficient_Sample_Analysis("Signal_4TeV_decayvx", Signal_4TeV_noPU, (0,35), "decayvx")
	#Efficient_Sample_Analysis("BG_decayvx", Background_noPU, (0,0.2), "decayvx")
	tfile1 = rt.TFile.Open("Thesis_Plots/root_files/Signal_2TeV_decayvx.root")
	decayvx_hist_2TeV = tfile1.Get("decay_vertex_R")
	tfile2 = rt.TFile.Open("Thesis_Plots/root_files/Signal_4TeV_decayvx.root")
	decayvx_hist_4TeV = tfile2.Get("decay_vertex_R")
	tfile3 = rt.TFile.Open("Thesis_Plots/root_files/BG_decayvx.root")
	decayvx_hist_BG = tfile3.Get("decay_vertex_R")
	#Histograms = [(decayvx_hist_2TeV,"2TeV-Signal"),(decayvx_hist_4TeV,"4TeV-Signal")]
	#DrawHistograms(Histograms, (0,35), "decay_vertex_R_Signal", "decay vertex R (cm)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	#DrawHistograms([(decayvx_hist_BG,"Background")], (0,0.2), "decay_vertex_R_Background", "decay vertex R (cm)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	
	decayvx_plots(decayvx_hist_2TeV,decayvx_hist_4TeV,decayvx_hist_BG)
	'''
	'''
	#Efficient_Sample_Analysis("Signal_2TeV_CSV", Signal_2TeV_noPU, (0,1), "CSV")
	#Efficient_Sample_Analysis("Signal_4TeV_CSV", Signal_4TeV_noPU, (0,1), "CSV")
	#Efficient_Sample_Analysis("BG_CSV", Background_noPU, (0,1), "CSV")
	tfile1 = rt.TFile.Open("Thesis_Plots/root_files/Signal_2TeV_CSV.root")
	CSV_hist_2TeV = tfile1.Get("CSV")
	tfile2 = rt.TFile.Open("Thesis_Plots/root_files/Signal_4TeV_CSV.root")
	CSV_hist_4TeV = tfile2.Get("CSV")
	tfile3 = rt.TFile.Open("Thesis_Plots/root_files/BG_CSV.root")
	CSV_hist_BG = tfile3.Get("CSV")
	Histograms = [(CSV_hist_BG,"background"),(CSV_hist_2TeV,"2TeV signal"),(CSV_hist_4TeV,"4TeV signal")]
	DrawHistograms(Histograms, (0,1), "CSV_b-tag", "p(b-jet)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''
	'''
	M2TeV_file = rt.TFile.Open("histogram_files/CheckInvariantMass_2TeV.root"); M2TeV = M2TeV_file.Get("M")
	M4TeV_file = rt.TFile.Open("histogram_files/CheckInvariantMass_4TeV.root"); M4TeV = M4TeV_file.Get("M")
	Histograms = [(M2TeV,"2TeV signal",3),(M4TeV,"4TeV signal",4)]
        DrawHistograms(Histograms, (0,5000), "Invariant_Mass", "M (GeV)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''


	#nice plot illustrating cluster matching

	'''	
	#with open("Grid.pkl",) as f:   #open coordinates of DetUnits for visual reference
        #       Grid = pickle.load(f)
	#ax = FCM.Initialize3DPlot('', 'x (cm)', 'y (cm)', 'z (cm)', grid=Grid)
	#HitClusters = CM.ClusterMatch(rt.TFile.Open(Signal_4TeV_noPU_String.format(20)), 0.1, 350, HadronsNotQuarks=True, Plot=True, Axes=ax, Save=False, dR_dist=False, LayerHist=False, EarlyBreak=5)
	CM.ClusterMatch2DPlot(rt.TFile.Open(Signal_4TeV_noPU_String.format(20)), 0.1, 350, HadronsNotQuarks=True, grid=np.load('Grid2D.npy'), EarlyBreak=5)
	'''

	#satisfactory example of decay after first layer falsely reconstructed by CSV	(obsolete)
	'''
	ax = FCM.Initialize3DPlot('SV Misalignment', 'x', 'y', 'z', grid=None)	
	X = [-1.5856226682662964, -1.5129914283752441, -1.4635052680969238, -1.463505744934082, -1.4635061025619507, -1.4563225507736206, -1.456322193145752, -1.3988516330718994, -3.379647731781006, -3.5370419025421143, -3.546525716781616, -3.5509400367736816, -3.553328514099121, -3.54559588432312, -3.558265209197998, -3.5627031326293945, -3.6817047595977783, -3.9403369426727295, -6.279918193817139, -8.275765419006348, -8.118766784667969, -8.110627174377441, -8.085073471069336, -8.053303718566895, -8.031609535217285, -8.011832237243652, -8.017459869384766, -7.996641159057617, -7.987757205963135, -7.982624053955078, -7.87465238571167, -7.837876796722412, -7.788309097290039, -7.677995681762695]
	Y = [-2.2767174243927, -2.347062826156616, -2.394984483718872, -2.394984722137451, -2.3949851989746094, -2.4019417762756348, -2.4019415378570557, -2.457594156265259, -6.061110496520996, -5.965805530548096, -5.960059642791748, -5.957389831542969, -5.955942153930664, -5.960626602172852, -5.952953338623047, -5.950264930725098, -5.878203392028809, -5.721562385559082, -9.087099075317383, -13.477447509765625, -13.571651458740234, -13.576526641845703, -13.591850280761719, -13.610928535461426, -13.62394905090332, -13.635809898376465, -13.632439613342285, -13.64491081237793, -13.650251388549805, -13.653338432312012, -13.718134880065918, -13.740201950073242, -13.76994800567627, -13.836077690124512]
	Z = [-10.43083381652832, -10.15747356414795, -10.228289604187012, -10.204048156738281, -10.181198120117188, -10.180895805358887, -10.202410697937012, -10.329291343688965, -15.388201713562012, -15.107477188110352, -15.201056480407715, -15.100199699401855, -15.138250350952148, -15.077401161193848, -15.110121726989746, -15.155170440673828, -15.075479507446289, -15.822460174560547, -21.128389358520508, -25.612314224243164, -25.637081146240234, -25.751169204711914, -25.867238998413086, -25.68181800842285, -25.633235931396484, -25.72307014465332, -25.644866943359375, -25.897769927978516, -25.779296875, -25.68834686279297, -25.583677291870117, -25.574443817138672, -25.540725708007812, -26.318405151367188]
	ax.scatter(X,Y,Z,color = 'blue',s=5)	
	MatchedTracks =  TM.TrackMatch(rt.TFile.Open(Signal_4TeV_noPU_String.format(20)), 0.1, 350, HadronsNotQuarks=True, Plot=True, Axes=ax, Save=False, dR_dist = False, EarlyBreak=40,manual_selection = (19,1056)	)		#20: (19,1056) best so far
	#changed the color of the trajectory tmporarily to green and shortened it to 0.8 of its length
	plt.show()
	'''

	#maybe plot showing detector modules in position space (a bit unnecessary, not?)

	
	#dR-discussion plots: dR-distribution and dR-dependng ROCs

	
	#hits per layer analysis: 'global' hits per layer, 'local' histograms
	'''
	Efficient_Layer_Hist2("4TeV",Signal_4TeV_noPU,0.1,minR=0,minPT1=200,minPT2=1000,Save=True)
	Efficient_Layer_Hist2("2TeV",Signal_2TeV_noPU,0.1,minR=0,minPT1=200,minPT2=1000,Save=True)
	Efficient_Layer_Hist2("Background",Background_noPU,0.1,minR=0,minPT1=200,minPT2=1000,Save=True)
	#FCM.Efficient_SeparateLayerHist([(Signal_2TeV_noPU,"2TeV"),(Signal_4TeV_noPU,"4TeV"),(Background_noPU,"Background")], (0,30), 0.1 , minPT=200,jet_pT=True, Save=True)
	'''
	
	#comparison of delta- and ratio-taggers for different layer combinations (ROC-Curves)
	'''
	General_Make_ROC_histograms('BG', Background_noPU,0.04, 200)	
	General_Make_ROC_histograms('2TeV', Signal_2TeV_noPU,0.04, 200)
	General_Make_ROC_histograms('4TeV', Signal_4TeV_noPU,0.04, 200)  	
	General_Make_ROC_histograms('BG', Background_noPU,0.06, 200)	
        General_Make_ROC_histograms('2TeV', Signal_2TeV_noPU,0.06, 200)
        General_Make_ROC_histograms('4TeV', Signal_4TeV_noPU,0.06, 200)
	General_Make_ROC_histograms('BG', Background_noPU,0.08, 200)	
        General_Make_ROC_histograms('2TeV', Signal_2TeV_noPU,0.08, 200)
        General_Make_ROC_histograms('4TeV', Signal_4TeV_noPU,0.08, 200)
	General_Make_ROC_histograms('BG', Background_noPU,0.1, 200)	
        General_Make_ROC_histograms('2TeV', Signal_2TeV_noPU,0.1, 200)
        General_Make_ROC_histograms('4TeV', Signal_4TeV_noPU,0.1, 200)
	General_Make_ROC_histograms('BG', Background_noPU,0.16, 200)	
        General_Make_ROC_histograms('2TeV', Signal_2TeV_noPU,0.16, 200)
        General_Make_ROC_histograms('4TeV', Signal_4TeV_noPU,0.16, 200)
	'''
	'''
	if len(sys.argv) > 1: dR = sys.argv[1] #most efficient way to make ROC curves for all values of dR is by putting it in externally as additional parameter
	tfile_4TeV= rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_dR{}.root".format(dR))
	tfile_BG = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_dR{}.root".format(dR))
	
	signal_4TeV_CSV = tfile_4TeV.Get("CSV")
	bg_CSV = tfile_BG.Get("CSV")
		
	signal_4TeV_ratio_21 = tfile_4TeV.Get("L2_L1")
	signal_4TeV_ratio_31 = tfile_4TeV.Get("L3_L1")
	signal_4TeV_ratio_41 = tfile_4TeV.Get("L4_L1")
	signal_4TeV_ratio_32 = tfile_4TeV.Get("L3_L2")
	signal_4TeV_ratio_42 = tfile_4TeV.Get("L4_L2")
	signal_4TeV_ratio_43 = tfile_4TeV.Get("L4_L3")
	bg_ratio_21 = tfile_BG.Get("L2_L1")
	bg_ratio_31 = tfile_BG.Get("L3_L1")
	bg_ratio_41 = tfile_BG.Get("L4_L1")
	bg_ratio_32 = tfile_BG.Get("L3_L2")
	bg_ratio_42 = tfile_BG.Get("L4_L2")
	bg_ratio_43 = tfile_BG.Get("L4_L3")
	signal_4TeV_diff_21 = tfile_4TeV.Get("L2-L1")
	signal_4TeV_diff_31 = tfile_4TeV.Get("L3-L1")
	signal_4TeV_diff_41 = tfile_4TeV.Get("L4-L1")
	signal_4TeV_diff_32 = tfile_4TeV.Get("L3-L2")
	signal_4TeV_diff_42 = tfile_4TeV.Get("L4-L2")
	signal_4TeV_diff_43 = tfile_4TeV.Get("L4-L3")
	bg_diff_21 = tfile_BG.Get("L2-L1")
	bg_diff_31 = tfile_BG.Get("L3-L1")
	bg_diff_41 = tfile_BG.Get("L4-L1")
	bg_diff_32 = tfile_BG.Get("L3-L2")
	bg_diff_42 = tfile_BG.Get("L4-L2")
	bg_diff_43 = tfile_BG.Get("L4-L3")

	bg_ZeroDiv_21, bg_ZeroDiv_31, bg_ZeroDiv_41, bg_ZeroDiv_32, bg_ZeroDiv_42, bg_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_BG_dR{}.csv".format(dR),delimiter=',')
	signal_4TeV_ZeroDiv_21, signal_4TeV_ZeroDiv_31, signal_4TeV_ZeroDiv_41, signal_4TeV_ZeroDiv_32, signal_4TeV_ZeroDiv_42, signal_4TeV_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_4TeV_dR{}.csv".format(dR),delimiter=',')

	ratio_histlist = [(signal_4TeV_ratio_21, bg_ratio_21, "L2/L1", "ratio", signal_4TeV_ZeroDiv_21, bg_ZeroDiv_21), (signal_4TeV_ratio_31, bg_ratio_31, "L3/L1", "ratio", signal_4TeV_ZeroDiv_31, bg_ZeroDiv_31), (signal_4TeV_ratio_41, bg_ratio_41, "L4/L1", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41), (signal_4TeV_ratio_32, bg_ratio_32, "L3/L2", "ratio", signal_4TeV_ZeroDiv_32, bg_ZeroDiv_32), (signal_4TeV_ratio_42, bg_ratio_42, "L4/L2", "ratio", signal_4TeV_ZeroDiv_42, bg_ZeroDiv_42), (signal_4TeV_ratio_43, bg_ratio_43, "L4/L3", "ratio", signal_4TeV_ZeroDiv_43, bg_ZeroDiv_43)]

	#ratio_histlist = [(signal_4TeV_ratio_21, bg_ratio_21, "L2/L1", "ratio", signal_4TeV_ZeroDiv_21, bg_ZeroDiv_21), (signal_4TeV_ratio_31, bg_ratio_31, "L3/L1", "ratio", signal_4TeV_ZeroDiv_31, bg_ZeroDiv_31), (signal_4TeV_ratio_41, bg_ratio_41, "L4/L1", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41), (signal_4TeV_ratio_32, bg_ratio_32, "L3/L2", "ratio", signal_4TeV_ZeroDiv_32, bg_ZeroDiv_32), (signal_4TeV_ratio_42, bg_ratio_42, "L4/L2", "ratio", signal_4TeV_ZeroDiv_42, bg_ZeroDiv_42), (signal_4TeV_ratio_43, bg_ratio_43, "L4/L3", "ratio", signal_4TeV_ZeroDiv_43, bg_ZeroDiv_43),(signal_4TeV_CSV, bg_CSV, "CSV", "CSV", 0, 0)]

	diff_histlist = [(signal_4TeV_diff_21, bg_diff_21, "L2-L1", "diff", 0, 0), (signal_4TeV_diff_31, bg_diff_31, "L3-L1", "diff", 0, 0), (signal_4TeV_diff_41, bg_diff_41, "L4-L1", "diff", 0, 0), (signal_4TeV_diff_32, bg_diff_32, "L3-L2", "diff", 0, 0), (signal_4TeV_diff_42, bg_diff_42, "L4-L2", "diff", 0, 0), (signal_4TeV_diff_43, bg_diff_43, "L4-L3", "diff", 0, 0)]

	#diff_histlist = [(signal_4TeV_diff_21, bg_diff_21, "L2-L1", "diff", 0, 0), (signal_4TeV_diff_31, bg_diff_31, "L3-L1", "diff", 0, 0), (signal_4TeV_diff_41, bg_diff_41, "L4-L1", "diff", 0, 0), (signal_4TeV_diff_32, bg_diff_32, "L3-L2", "diff", 0, 0), (signal_4TeV_diff_42, bg_diff_42, "L4-L2", "diff", 0, 0), (signal_4TeV_diff_43, bg_diff_43, "L4-L3", "diff", 0, 0),(signal_4TeV_CSV, bg_CSV, "CSV", "CSV", 0, 0)]

	dR_string = str(dR).translate(None,'.')
	General_Make_ROC_Curves("ratio_taggers_dR{}".format(dR_string), ratio_histlist,log=True,dR=dR)
	General_Make_ROC_Curves("diff_taggers_dR{}".format(dR_string), diff_histlist,log=True,dR=dR)
	'''
	#Discriminant Histograms corresponding to all combinatinos of discriminants: ratios at dR<0.1 and deltas at dR<0.04
	'''
	tfile_2TeV_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_2TeV_dR{}.root".format(0.04))
	tfile_4TeV_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_dR{}.root".format(0.04))
	tfile_BG_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_dR{}.root".format(0.04))
	tfile_2TeV_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_2TeV_dR{}.root".format(0.1))
	tfile_4TeV_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_dR{}.root".format(0.1))
	tfile_BG_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_dR{}.root".format(0.1))
		
	ratio_L2_L1_2TeV = tfile_2TeV_01.Get("L2_L1"); ratio_L3_L1_2TeV = tfile_2TeV_01.Get("L3_L1"); ratio_L4_L1_2TeV = tfile_2TeV_01.Get("L4_L1")
	ratio_L3_L2_2TeV = tfile_2TeV_01.Get("L3_L2"); ratio_L4_L2_2TeV = tfile_2TeV_01.Get("L4_L2"); ratio_L4_L3_2TeV = tfile_2TeV_01.Get("L4_L3")
	ratio_L2_L1_4TeV = tfile_4TeV_01.Get("L2_L1"); ratio_L3_L1_4TeV = tfile_4TeV_01.Get("L3_L1"); ratio_L4_L1_4TeV = tfile_4TeV_01.Get("L4_L1")
	ratio_L3_L2_4TeV = tfile_4TeV_01.Get("L3_L2"); ratio_L4_L2_4TeV = tfile_4TeV_01.Get("L4_L2"); ratio_L4_L3_4TeV = tfile_4TeV_01.Get("L4_L3")
	ratio_L2_L1_BG = tfile_BG_01.Get("L2_L1"); ratio_L3_L1_BG = tfile_BG_01.Get("L3_L1"); ratio_L4_L1_BG = tfile_BG_01.Get("L4_L1")
	ratio_L3_L2_BG = tfile_BG_01.Get("L3_L2"); ratio_L4_L2_BG = tfile_BG_01.Get("L4_L2"); ratio_L4_L3_BG = tfile_BG_01.Get("L4_L3")
	diff_L2_L1_2TeV = tfile_2TeV_01.Get("L2-L1"); diff_L3_L1_2TeV = tfile_2TeV_01.Get("L3-L1"); diff_L4_L1_2TeV = tfile_2TeV_01.Get("L4-L1")
	diff_L3_L2_2TeV = tfile_2TeV_01.Get("L3-L2"); diff_L4_L2_2TeV = tfile_2TeV_01.Get("L4-L2"); diff_L4_L3_2TeV = tfile_2TeV_01.Get("L4-L3")
	diff_L2_L1_4TeV = tfile_4TeV_01.Get("L2-L1"); diff_L3_L1_4TeV = tfile_4TeV_01.Get("L3-L1"); diff_L4_L1_4TeV = tfile_4TeV_01.Get("L4-L1")
	diff_L3_L2_4TeV = tfile_4TeV_01.Get("L3-L2"); diff_L4_L2_4TeV = tfile_4TeV_01.Get("L4-L2"); diff_L4_L3_4TeV = tfile_4TeV_01.Get("L4-L3")
	diff_L2_L1_BG = tfile_BG_01.Get("L2-L1"); diff_L3_L1_BG = tfile_BG_01.Get("L3-L1"); diff_L4_L1_BG = tfile_BG_01.Get("L4-L1")
	diff_L3_L2_BG = tfile_BG_01.Get("L3-L2"); diff_L4_L2_BG = tfile_BG_01.Get("L4-L2"); diff_L4_L3_BG = tfile_BG_01.Get("L4-L3")

	DrawHistograms([(ratio_L2_L1_BG,'background'),(ratio_L2_L1_2TeV,'2TeV signal'),(ratio_L2_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio1', 'L2/L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(ratio_L3_L1_BG,'background'),(ratio_L3_L1_2TeV,'2TeV signal'),(ratio_L3_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio2', 'L3/L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(ratio_L4_L1_BG,'background'),(ratio_L4_L1_2TeV,'2TeV signal'),(ratio_L4_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio3', 'L4/L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(ratio_L3_L2_BG,'background'),(ratio_L3_L2_2TeV,'2TeV signal'),(ratio_L3_L2_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio4', 'L3/L2', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(ratio_L4_L2_BG,'background'),(ratio_L4_L2_2TeV,'2TeV signal'),(ratio_L4_L2_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio5', 'L4/L2', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(ratio_L4_L3_BG,'background'),(ratio_L4_L3_2TeV,'2TeV signal'),(ratio_L4_L3_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_ratio6', 'L4/L3', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L2_L1_BG,'background'),(diff_L2_L1_2TeV,'2TeV signal'),(diff_L2_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta1', 'L2-L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L3_L1_BG,'background'),(diff_L3_L1_2TeV,'2TeV signal'),(diff_L3_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta2', 'L3-L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L4_L1_BG,'background'),(diff_L4_L1_2TeV,'2TeV signal'),(diff_L4_L1_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta3', 'L4-L1', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L3_L2_BG,'background'),(diff_L3_L2_2TeV,'2TeV signal'),(diff_L3_L2_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta4', 'L3-L2', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L4_L2_BG,'background'),(diff_L4_L2_2TeV,'2TeV signal'),(diff_L4_L2_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta5', 'L4-L2', '(a.u)', Save=True, Normalize=True)
	DrawHistograms([(diff_L4_L3_BG,'background'),(diff_L4_L3_2TeV,'2TeV signal'),(diff_L4_L3_4TeV,'4TeV signal')], (0,10), 'all_Discriminant_hists/Discriminant_Hist_delta6', 'L4-L3', '(a.u)', Save=True, Normalize=True)
	'''
	
	#best ROC-curves global and single pT threshold (discriminant variables from above)
	
	'''	
	bg_ZeroDiv_21, bg_ZeroDiv_31, bg_ZeroDiv_41, bg_ZeroDiv_32, bg_ZeroDiv_42, bg_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_BG_dR{}.csv".format(0.1),delimiter=',')
	signal_4TeV_ZeroDiv_21, signal_4TeV_ZeroDiv_31, signal_4TeV_ZeroDiv_41, signal_4TeV_ZeroDiv_32, signal_4TeV_ZeroDiv_42, signal_4TeV_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_4TeV_dR{}.csv".format(0.1),delimiter=',')
	
	
	best_histlist = [(signal_4TeV_ratio, bg_ratio, "L4/L1", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41), (signal_4TeV_diff, bg_diff, "L4-L1", "diff", 0,0),(signal_4TeV_CSV, bg_CSV, "CSV", "CSV", 0,0)]
	#General_Make_ROC_Curves("Best_Taggers", best_histlist,log=True,print_cut=True)
	
	#General_Make_ROC_histograms('BG_HPT', Background_noPU,0.04, 1200)	
	#General_Make_ROC_histograms('4TeV_HPT', Signal_4TeV_noPU,0.04, 1200)  	
	#General_Make_ROC_histograms('BG_HPT', Background_noPU,0.1, 1200)	
        #General_Make_ROC_histograms('4TeV_HPT', Signal_4TeV_noPU,0.1, 1200)
	
	tfile_4TeV_004_HPT = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_HPT_dR{}.root".format(0.04))
	tfile_BG_004_HPT = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_HPT_dR{}.root".format(0.04))
	tfile_4TeV_01_HPT = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_HPT_dR{}.root".format(0.1))
	tfile_BG_01_HPT = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_HPT_dR{}.root".format(0.1))
	signal_4TeV_ratio_HPT = tfile_4TeV_01_HPT.Get("L4_L1")
	bg_ratio_HPT = tfile_BG_01_HPT.Get("L4_L1")
	signal_4TeV_diff_HPT = tfile_4TeV_004_HPT.Get("L4-L1")
	bg_diff_HPT = tfile_BG_004_HPT.Get("L4-L1")
	signal_4TeV_CSV_HPT = tfile_4TeV_01_HPT.Get("CSV")	
	bg_CSV_HPT = tfile_BG_01_HPT.Get("CSV")	

	bg_ZeroDiv_21_HPT, bg_ZeroDiv_31_HPT, bg_ZeroDiv_41_HPT, bg_ZeroDiv_32_HPT, bg_ZeroDiv_42_HPT, bg_ZeroDiv_43_HPT = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_BG_HPT_dR{}.csv".format(0.1),delimiter=',')
	signal_4TeV_ZeroDiv_21_HPT, signal_4TeV_ZeroDiv_31_HPT, signal_4TeV_ZeroDiv_41_HPT, signal_4TeV_ZeroDiv_32_HPT, signal_4TeV_ZeroDiv_42_HPT, signal_4TeV_ZeroDiv_43_HPT = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_4TeV_HPT_dR{}.csv".format(0.1),delimiter=',')

	diff_HPT_histlist = [(signal_4TeV_diff, bg_diff, "L4-L1 ($p_T$>200GeV)", "diff", 0,0, 'green','-'),(signal_4TeV_diff_HPT, bg_diff_HPT, "L4-L1 ($p_T$>1200GeV)", "diff", 0,0, 'green','--'),(signal_4TeV_CSV, bg_CSV, "CSV ($p_T$>200GeV)", "CSV", 0,0,'blue','-'),(signal_4TeV_CSV_HPT, bg_CSV_HPT, "CSV ($p_T$>1200GeV)", "CSV", 0,0,'blue','--')]
	ratio_HPT_histlist = [(signal_4TeV_ratio, bg_ratio, "L4/L1 ($p_T$>200GeV)", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41, 'red','-'),(signal_4TeV_ratio_HPT, bg_ratio_HPT, "L4/L1 ($p_T$>1200GeV)", "ratio", signal_4TeV_ZeroDiv_41_HPT, bg_ZeroDiv_41_HPT, 'red', '--'),(signal_4TeV_CSV, bg_CSV, "CSV ($p_T$>200GeV)", "CSV", 0,0, 'blue','-'),(signal_4TeV_CSV_HPT, bg_CSV_HPT, "CSV ($p_T$>1200GeV)", "CSV", 0,0, 'blue', '--')]
	
	#General_Make_ROC_Curves("diff_Taggers_high_pT", diff_HPT_histlist,log=True,print_cut=True,dR=0.04)
	#General_Make_ROC_Curves("ratio_Taggers_high_pT", ratio_HPT_histlist,log=True,print_cut=True,dR=0.1)
	'''
	
	#Discriminant histograms above pT threshold
	'''	
	DrawHistograms([(bg_ratio_HPT,'background'),(signal_4TeV_ratio_HPT,'4TeV signal',4)], (0,10), 'L4_L1_Discriminant_Hist_high_pT', 'L4/L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	DrawHistograms([(bg_diff_HPT,'background'),(signal_4TeV_diff_HPT,'4TeV signal',4)], (-22,22), 'L4-L1_Discriminant_Hist_high_pT', 'L4-L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	DrawHistograms([(bg_CSV_HPT,'background'),(signal_4TeV_CSV_HPT,'4TeV signal',4)], (0,1), 'CSV_Discriminant_Hist_high_pT', 'p(b-jet)', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	'''
	
	#efficiency vs cut (histograms and ZeroDiv form above)
	'''
	efficiency_vs_cut("L4_L1", signal_4TeV_ratio, bg_ratio, signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41, (0,6), 38)
	efficiency_vs_cut("L4-L1", signal_4TeV_diff, bg_diff, 0, 0, (-22,22), 43)
	efficiency_vs_cut("CSV", signal_4TeV_CSV, bg_CSV, 0, 0, (0,1), 59)
	'''

	#tagged jets vs pT (also exclusive)
	'''
	efficient_tagged_jets_hist([(Signal_4TeV_noPU,"4TeV_signal",(0,2500))],"L4_L1", 1.833, 0.66667, 60, Difference=False, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal_4TeV_noPU,"4TeV_signal",(0,2500))],"L4-L1", 5, 0.66667, 60, Difference=True, mode="pT_jet",Save=True)
	efficient_tagged_jets_hist([(Signal_4TeV_noPU,"4TeV_signal",(0,30))],"L4_L1", 1.833, 0.66667, 60, Difference=False, mode="decay_vx",Save=True)
	efficient_tagged_jets_hist([(Signal_4TeV_noPU,"4TeV_signal",(0,30))],"L4-L1", 5, 0.66667, 60, Difference=True, mode="decay_vx",Save=True)
	'''
	'''
	tagged_diff_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_pT_jet_4TeV_signalL4-L1.root')
	tagged_ratio_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_pT_jet_4TeV_signalL4_L1.root') #just switch between source file to choose between pT or decayvx
	#tagged_diff_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_decay_vx_4TeV_signalL4-L1.root')
	#tagged_ratio_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_decay_vx_4TeV_signalL4_L1.root')
	AllJets_thist = tagged_diff_file.Get('4TeV_signal_AllJets')
	CSV_thist = tagged_diff_file.Get('4TeV_signal_CSV')
	diff_thist = tagged_diff_file.Get('4TeV_signal_Discriminant')
	ratio_thist = tagged_ratio_file.Get('4TeV_signal_Discriminant')
	
	#DrawHistograms([(AllJets_thist,"all jets"), (ratio_thist, "L4/L1"), (diff_thist, "L4-L1"), (CSV_thist,"CSV")], (0,2500), "tagged_jets_vs_pT", 'jet-pT', "# jets", Save=True,Normalize=False)
	#DrawHistograms([(AllJets_thist,"all jets"), (ratio_thist, "L4/L1"), (diff_thist, "L4-L1"), (CSV_thist,"CSV")], (0,30), "tagged_jets_vs_decayvx", 'decay vertex R', "# jets", Save=True,Normalize=False)


#	binned_exclusive_tagged_jets_hist("L4_L1_single_cut",Signal_4TeV_noPU, "L4/L1", [1.833], [0.66667], [0,2500], (0,2500), 60, Difference=False, mode="pT_jet", y_max = 1300, Save=True, Stacked=True, AllJets=AllJets_thist)
	#binned_exclusive_tagged_jets_hist("L4-L1_single_cut",Signal_4TeV_noPU, "L4-L1", [5], [0.66667], [0,2500], (0,2500), 60, Difference=True, mode="pT_jet", y_max = 1300, Save=True, Stacked=True, AllJets=AllJets_thist)
	'''

	#efficiency vs pT
	
	plot_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500] #used for a nice looking plot
	'''
	#efficient_tagged_jets_hist([(Background_noPU,"BG",(0,2500))],"L4_L1", 1.833, 0.66667, 60, Difference=False, mode="pT_jet",Save=True)
	#efficient_tagged_jets_hist([(Background_noPU,"BG",(0,2500))],"L4-L1", 5, 0.66667, 60, Difference=True, mode="pT_jet",Save=True)
	bg_tagged_diff_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_pT_jet_BGL4-L1.root')
        bg_tagged_ratio_file = rt.TFile.Open('Thesis_Plots/root_files/tagged_jets_vs_pT_jet_BGL4_L1.root')
	bg_AllJets_thist = 	FCM.RebinHist(bg_tagged_diff_file.Get('BG_AllJets'),"AllJets",plot_bins)
        bg_CSV_thist = 		FCM.RebinHist(bg_tagged_diff_file.Get('BG_CSV'),"CSV",plot_bins)
        bg_diff_thist = 	FCM.RebinHist(bg_tagged_diff_file.Get('BG_Discriminant'),"L4-L1",plot_bins)
        bg_ratio_thist = 	FCM.RebinHist(bg_tagged_ratio_file.Get('BG_Discriminant'),"L4/L1",plot_bins)

	AllJets_thist2 = 	FCM.RebinHist(AllJets_thist,"AllJets",plot_bins)
        CSV_thist2 = 		FCM.RebinHist(CSV_thist,"CSV",plot_bins)
        diff_thist2 = 		FCM.RebinHist(diff_thist,"L4-L1",plot_bins)
        ratio_thist2 = 		FCM.RebinHist(ratio_thist,"L4/L1",plot_bins)

	Efficiency_vs_pT("4TeV-signal",[(ratio_thist2,"L4/L1"),(diff_thist2,"L4-L1"),(CSV_thist2,"CSV")], AllJets_thist2, 0.6, Save=True,legend_shift=True)
	Efficiency_vs_pT("Background",[(bg_ratio_thist,"L4/L1"),(bg_diff_thist,"L4-L1"),(bg_CSV_thist,"CSV")], bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True)
	'''
	'''
	exclusive_tagged_jets_hist("4TeV", Signal_4TeV_noPU , "L4/L1", 1.833, 0.667,(200,2500), 60, Difference=False, mode="pT_jet",Save=True)
	exclusive_tagged_jets_hist("4TeV", Signal_4TeV_noPU , "L4-L1", 5, 0.667,(200,2500), 60, Difference=True, mode="pT_jet",Save=True)
	'''

	#study on pT-bins

	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500] #used for the pT-depending cuts
	coarse_bins = [0,1200,1800,2500]	
	
	#efficient_Make_Binned_ROC_histograms('Signal', Signal_both_noPU, cut_bins)    
        #efficient_Make_Binned_ROC_histograms('BG', Background_noPU, cut_bins)
	
	#FCM.find_cuts('Thesis_Plots/root_files/BG_histograms.root', cut_bins, ZeroDiv_path="Thesis_Plots/root_files/BG_ZeroDiv.csv")

	delta_cuts_noPU = [3, 5, 6, 6, 7, 8, 8, 9, 9]
        ratio_cuts_noPU = [1.833, 1.833, 2.0, 2.0, 2.0, 2.167, 2.167, 2.167, 2.167]
        CSV_cuts_noPU = [0.633, 0.65, 0.667, 0.683, 0.7, 0.717, 0.733, 0.767, 0.767]
	'''
	#efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,2500))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        #efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,2500))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	#efficient_binned_tagged_jets_hist([(Background_noPU, "BG",(0,2500))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        #efficient_binned_tagged_jets_hist([(Background_noPU, "BG",(0,2500))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	
	
	
	binned_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_SignalL4-L1.root")
        binned_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_SignalL4_L1.root")
	binned_bg_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BGL4-L1.root")
        binned_bg_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BGL4_L1.root")
	binned_AllJets_thist = 		FCM.RebinHist(binned_tagged_diff_file.Get('Signal_AllJets'),"AllJets",plot_bins)			#binned_tagged_diff_file.Get('Signal_AllJets')		
        binned_CSV_thist = 		FCM.RebinHist(binned_tagged_diff_file.Get('Signal_CSV'),"CSV",plot_bins)                       	#binned_tagged_diff_file.Get('Signal_CSV')		
        binned_diff_thist = 		FCM.RebinHist(binned_tagged_diff_file.Get('Signal_Discriminant'),"L4-L1",plot_bins)            	#binned_tagged_diff_file.Get('Signal_Discriminant')	
        binned_ratio_thist = 		FCM.RebinHist(binned_tagged_ratio_file.Get('Signal_Discriminant'),"L4/L1",plot_bins)           	#binned_tagged_ratio_file.Get('Signal_Discriminant')	
	binned_bg_AllJets_thist = 	FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_AllJets'),"AllJets",plot_bins)                	#binned_bg_tagged_diff_file.Get('BG_AllJets')		
        binned_bg_CSV_thist = 		FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_CSV'),"CSV",plot_bins)                        	#binned_bg_tagged_diff_file.Get('BG_CSV')		
        binned_bg_diff_thist = 		FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_Discriminant'),"L4-L1",plot_bins)             	#binned_bg_tagged_diff_file.Get('BG_Discriminant')	
        binned_bg_ratio_thist = 	FCM.RebinHist(binned_bg_tagged_ratio_file.Get('BG_Discriminant'),"L4/L1",plot_bins)            	#binned_bg_tagged_ratio_file.Get('BG_Discriminant')	
	

	#binned_exclusive_tagged_jets_hist("Signal_binned",Signal_both_noPU, "L4/L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, (200,2500), 60, Difference=False, mode="pT_jet",Save=True, Stacked=True, AllJets = binned_AllJets_thist)
	#binned_exclusive_tagged_jets_hist("Signal_binned",Signal_both_noPU, "L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, (200,2500), 60, Difference=True, mode="pT_jet",Save=True, Stacked=True, AllJets = binned_AllJets_thist)
	
	#Relative_Gain_Plots_binned()
	
	Efficiency_vs_pT("Signal_binned",[(binned_ratio_thist,"L4/L1"),(binned_diff_thist,"L4-L1"),(binned_CSV_thist,"CSV")], binned_AllJets_thist, 0.6, Save=True,legend_shift=True)
        Efficiency_vs_pT("Background_binned",[(binned_bg_ratio_thist,"L4/L1"),(binned_bg_diff_thist,"L4-L1"),(binned_bg_CSV_thist,"CSV")], binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True)
	
	#DrawHistograms([(binned_AllJets_thist,"all jets"), (binned_ratio_thist, "L4/L1"), (binned_diff_thist, "L4-L1"), (binned_CSV_thist,"CSV")], (0,2500), "tagged_jets_vs_pT_binned", 'jet-pT', "# jets", Save=True,Normalize=False)
	
	#efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,30))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="decay_vx",Save=True)
        #efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,30))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="decay_vx",Save=True)
	
	binned_dvx_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_decay_vx_SignalL4-L1.root")
        binned_dvx_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_decay_vx_SignalL4_L1.root")
	binned_dvx_AllJets_thist = 		FCM.RebinHist(binned_tagged_diff_file.Get('Signal_AllJets'),"AllJets",plot_bins)		#binned_dvx_tagged_diff_file.Get('Signal_AllJets')	
        binned_dvx_CSV_thist = 			FCM.RebinHist(binned_tagged_diff_file.Get('Signal_CSV'),"CSV",plot_bins)                       #binned_dvx_tagged_diff_file.Get('Signal_CSV')		
        binned_dvx_diff_thist = 		FCM.RebinHist(binned_tagged_diff_file.Get('Signal_Discriminant'),"L4-L1",plot_bins)            #binned_dvx_tagged_diff_file.Get('Signal_Discriminant')	
        binned_dvx_ratio_thist = 		FCM.RebinHist(binned_tagged_ratio_file.Get('Signal_Discriminant'),"L4/L1",plot_bins)           #binned_dvx_tagged_ratio_file.Get('Signal_Discriminant')	

	#DrawHistograms([(binned_dvx_AllJets_thist,"all jets"), (binned_dvx_ratio_thist, "L4/L1"), (binned_dvx_diff_thist, "L4-L1"), (binned_dvx_CSV_thist,"CSV")], (0,30), "tagged_jets_vs_decayvx_binned", 'decay vertex R', "# jets", Save=True,Normalize=False)
	
	
	#efficient_Make_Binned_ROC_histograms('Signal_coarse-binned', Signal_both_noPU, coarse_bins)
        #efficient_Make_Binned_ROC_histograms('BG_coarse-binned', Background_noPU, coarse_bins)
	#Make_Binned_ROC_Curves('coarse_binned_ratio','Signal_coarse-binned','BG_coarse-binned',coarse_bins, diff=False,log=True, dR=0.1)
	#Make_Binned_ROC_Curves('coarse_binned_diff','Signal_coarse-binned','BG_coarse-binned',coarse_bins, diff=True,log=True, dR=0.04)
	'''	
	
	#PU study of cut-based taggers
	'''
	#PU_histograms(Signal_2TeV_PU, Signal_4TeV_PU, Background_PU)	
	
	PU_file = rt.TFile.Open("Thesis_Plots/root_files/PU_distributions.root")
	PU_2TeV = PU_file.Get("PU_2TeV"); PU_4TeV = PU_file.Get("PU_4TeV"); PU_BG = PU_file.Get("PU_BG")
	DrawHistograms([(PU_BG,"background"), (PU_2TeV,"2TeV signal"), (PU_4TeV,"4TeV signal")], (0,80), "PU_distributions", '#PV', "(a.u.)", Save=True,Normalize=True)
	'''
	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500] #used for the pT-depending cuts
	'''
	efficient_Make_Binned_ROC_histograms('Signal_PU', Signal_both_PU, cut_bins)    
        efficient_Make_Binned_ROC_histograms('BG_PU', Background_PU, cut_bins)

	FCM.find_cuts('Thesis_Plots/root_files/BG_PU_histograms.root',cut_bins, ZeroDiv_path = "Thesis_Plots/root_files/BG_PU_ZeroDiv.csv")
	'''
	delta_cuts_PU = [4, 5, 6, 7, 8, 8, 9, 9, 11]
        ratio_cuts_PU = [1.833, 2.0, 2.0, 2.167, 2.167, 2.167, 2.167, 2.167, 2.167]
        CSV_cuts_PU = [0.65, 0.667, 0.683, 0.7, 0.717, 0.733, 0.75, 0.783, 0.783]	
	'''
	efficient_binned_tagged_jets_hist([(Signal_both_PU, "Signal_PU",(0,2500))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        efficient_binned_tagged_jets_hist([(Signal_both_PU, "Signal_PU",(0,2500))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Background_PU, "BG_PU",(0,2500))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        efficient_binned_tagged_jets_hist([(Background_PU, "BG_PU",(0,2500))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	'''
	'''
	PU_binned_tagged_diff_file = 		rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PUL4-L1.root")
        PU_binned_tagged_ratio_file = 		rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PUL4_L1.root")
	PU_binned_bg_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PUL4-L1.root")
        PU_binned_bg_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PUL4_L1.root")
	PU_binned_AllJets_thist = 		FCM.RebinHist(PU_binned_tagged_diff_file.Get('Signal_PU_AllJets'),"AllJets",plot_bins)		#PU_binned_tagged_diff_file.Get('Signal_AllJets')	
        PU_binned_CSV_thist = 			FCM.RebinHist(PU_binned_tagged_diff_file.Get('Signal_PU_CSV'),"CSV",plot_bins)			#PU_binned_tagged_diff_file.Get('Signal_CSV')		
        PU_binned_diff_thist = 			FCM.RebinHist(PU_binned_tagged_diff_file.Get('Signal_PU_Discriminant'),"L4-L1",plot_bins)	#PU_binned_tagged_diff_file.Get('Signal_Discriminant')	
        PU_binned_ratio_thist = 		FCM.RebinHist(PU_binned_tagged_ratio_file.Get('Signal_PU_Discriminant'),"L4/L1",plot_bins)	#PU_binned_tagged_ratio_file.Get('Signal_Discriminant')	
	PU_binned_bg_AllJets_thist = 		FCM.RebinHist(PU_binned_bg_tagged_diff_file.Get('BG_PU_AllJets'),"AllJets",plot_bins)		#PU_binned_bg_tagged_diff_file.Get('BG_AllJets')		
        PU_binned_bg_CSV_thist = 		FCM.RebinHist(PU_binned_bg_tagged_diff_file.Get('BG_PU_CSV'),"CSV",plot_bins)			#PU_binned_bg_tagged_diff_file.Get('BG_CSV')		
        PU_binned_bg_diff_thist = 		FCM.RebinHist(PU_binned_bg_tagged_diff_file.Get('BG_PU_Discriminant'),"L4-L1",plot_bins)	#PU_binned_bg_tagged_diff_file.Get('BG_Discriminant')	
        PU_binned_bg_ratio_thist = 		FCM.RebinHist(PU_binned_bg_tagged_ratio_file.Get('BG_PU_Discriminant'),"L4/L1",plot_bins)	#PU_binned_bg_tagged_ratio_file.Get('BG_Discriminant')	
	
	Efficiency_vs_pT("Signal_binned_PU",[(PU_binned_ratio_thist,"L4/L1"),(PU_binned_diff_thist,"L4-L1"),(PU_binned_CSV_thist,"CSV")], PU_binned_AllJets_thist, 0.65, Save=True,legend_shift=True)
        Efficiency_vs_pT("Background_binned_PU",[(PU_binned_bg_ratio_thist,"L4/L1"),(PU_binned_bg_diff_thist,"L4-L1"),(PU_binned_bg_CSV_thist,"CSV")], PU_binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True)
	'''
        '''
	binned_exclusive_tagged_jets_hist("Signal_binned_PU",Signal_both_PU, "L4/L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, (200,2500), 60, Difference=False, mode="pT_jet",Save=True)
        binned_exclusive_tagged_jets_hist("Signal_binned_PU",Signal_both_PU, "L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, (200,2500), 60, Difference=True, mode="pT_jet",Save=True)
        '''
	'''
	binned_efficiency_vs_PU('Efficiency', Signal_both_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.8, pT_Cut=200, pT_Mode="jet")
        binned_efficiency_vs_PU('Efficiency', Signal_both_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.8, pT_Cut=1200, pT_Mode="jet")

        binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.3, pT_Cut=200, pT_Mode="jet",BG=True)
        binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.3, pT_Cut=1200, pT_Mode="jet",BG=True)
	'''

	#comparison of different ANN models
	
	from sklearn.metrics import roc_curve
	import keras as kr
	print "packages imported"
	
	'''
	#noPU
	model_noPT = kr.models.load_model("Submitted_Models/model_noPU_functional.h5")
	model_withPT = kr.models.load_model("Submitted_Models/model_noPU_functional_withPT.h5")
	model_type = "functional"
	print "model loaded"

	x_data = np.load("Submitted_Models/data/noPU_both_withPT/test_x.npy")
	CSV = np.load("Submitted_Models/data/noPU_both_withPT/test_CSV.npy")
	labels = np.load("Submitted_Models/data/noPU_both_withPT/test_y.npy")
	pT = np.load("Submitted_Models/data/noPU_both_withPT/test_feature.npy")
	print "data loaded"
	'''	
	
	#withPU both
	model_noPT = kr.models.load_model("Submitted_Models/model_withPU_functional.h5")
	model_withPT = kr.models.load_model("Submitted_Models/model_withPU_functional_withPT.h5")
	model_withPV = kr.models.load_model("Submitted_Models/model_withPU_functional_withPV.h5")
        model_type = "functional"
	print "model loaded"
	
	x_data = np.load("Submitted_Models/data/withPU_both_withPT/test_x.npy")
	CSV = np.load("Submitted_Models/data/withPU_both_withPT/test_CSV.npy")
	labels = np.load("Submitted_Models/data/withPU_both_withPT/test_y.npy")
	pT = np.load("Submitted_Models/data/withPU_both_withPT/test_feature.npy")
	print "data loaded"
	
	'''	
	pred_y = model.predict(ANN_functional_shape(x_data))
	fpr, tpr, thr = roc_curve(np.array(labels),pred_y)
	fpr2, tpr2, thr2 = roc_curve(np.array(labels),CSV)
	plt.figure()
	plt.plot(tpr,1-fpr,label='ANN')
	plt.plot(tpr2,1-fpr2,label='CSV')
	plt.legend(loc=3)
	plt.show()
	'''
	
	signal_x_data = x_data[labels==1]
	signal_pT = pT[labels==1]
	signal_CSV = CSV[labels==1]
	bg_x_data = x_data[labels==0]
	bg_pT = pT[labels==0]
	bg_CSV = CSV[labels==0]
		
	
	#noPU
	'''
	#print "noPT:"
	#ANN_Make_Binned_ROC_histograms("Signal_ANN",model_noPT, signal_x_data, signal_pT, signal_CSV, cut_bins)
	#ANN_Make_Binned_ROC_histograms("BG_ANN",model_noPT, bg_x_data, bg_pT, bg_CSV, cut_bins)
	#FCM.find_cuts('Thesis_Plots/root_files/BG_ANN_histograms.root',cut_bins,ANN=True)
	#print "withPT:"
	#ANN_Make_Binned_ROC_histograms("Signal_ANN",model_withPT, signal_x_data, signal_pT, signal_CSV, cut_bins,addFeature="pT")
	#ANN_Make_Binned_ROC_histograms("BG_ANN",model_withPT, bg_x_data, bg_pT, bg_CSV, cut_bins,addFeature="pT")
	#FCM.find_cuts('Thesis_Plots/root_files/BG_ANN_histograms.root',cut_bins,ANN=True)
	
	#ANN_cuts = [0.55, 0.583, 0.6, 0.667, 0.633, 0.6, 0.7, 0.667, 0.567] #from model before (complex functional)
	ANN_noPT_cuts = [0.55, 0.567, 0.6, 0.633, 0.6, 0.617, 0.683, 0.6, 0.533]
	ANN_withPT_cuts = [0.633, 0.633, 0.633, 0.633, 0.617, 0.6, 0.6, 0.65, 0.55]
	CSV_cuts = [0.633, 0.65, 0.667, 0.7, 0.667, 0.7, 0.817, 0.817, 0.767]
	
	#datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_noPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_noPT", (0,2500))]
	#ANN_binned_tagged_jets_hist(datalist, model_noPT, ANN_noPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True)
	#datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_withPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_withPT", (0,2500))]
	#ANN_binned_tagged_jets_hist(datalist, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True,addFeature="pT")
	
	
	tagged_ANN_noPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_noPTANN.root")
	bg_tagged_ANN_noPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_noPTANN.root")
	binned_AllJets_thist = 		FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_noPT_AllJets'),"AllJets",plot_bins)   					#tagged_ANN_noPT_file.Get('Signal_noPT_AllJets')			
        binned_CSV_thist = 		FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_noPT_CSV'),"CSV",plot_bins)           					#tagged_ANN_noPT_file.Get('Signal_noPT_CSV')				
        binned_ANN_noPT_thist = 	FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_noPT_Discriminant'),"ANN",plot_bins)					#tagged_ANN_noPT_file.Get('Signal_noPT_Discriminant')			
	binned_bg_AllJets_thist = 	FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_noPT_AllJets'),"AllJets",plot_bins)    					#bg_tagged_ANN_noPT_file.Get('BG_noPT_AllJets')				
        binned_bg_CSV_thist = 		FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_noPT_CSV'),"CSV",plot_bins)            					#bg_tagged_ANN_noPT_file.Get('BG_noPT_CSV')				
        binned_bg_ANN_noPT_thist =	FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_noPT_Discriminant'),"ANN",plot_bins) 					#bg_tagged_ANN_noPT_file.Get('BG_noPT_Discriminant')			
	tagged_ANN_withPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_withPTANN.root")
	bg_tagged_ANN_withPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_withPTANN.root")
	binned_ANN_withPT_thist = 	FCM.RebinHist(tagged_ANN_withPT_file.Get('Signal_withPT_Discriminant'),"ANN",plot_bins)			#tagged_ANN_withPT_file.Get('Signal_withPT_Discriminant')				
	binned_bg_ANN_withPT_thist = 	FCM.RebinHist(bg_tagged_ANN_withPT_file.Get('BG_withPT_Discriminant'),"ANN",plot_bins) 			#bg_tagged_ANN_withPT_file.Get('BG_withPT_Discriminant')				
	
	#ANN_exclusive_tagged_jets_hist('Signal_ANN_noPT', model_noPT, signal_x_data, signal_pT, signal_CSV, ANN_noPT_cuts, CSV_cuts, cut_bins, (0,2500), 60, mode="pT_jet", y_max=870, Save=True, Stacked=True, AllJets=binned_AllJets_thist )
	#ANN_exclusive_tagged_jets_hist('Signal_ANN_withPT', model_withPT, signal_x_data, signal_pT, signal_CSV, ANN_withPT_cuts, CSV_cuts, cut_bins, (0,2500), 60, mode="pT_jet", y_max=870, Save=True,addFeature="pT", Stacked=True, AllJets=binned_AllJets_thist )

	coarse_bins = [0,1200,1500,2500]	
	#coarse_bins = [0,800,1000,1200]	
	#coarse_bins = [0,2500]
	#ANN_Make_Binned_ROC_histograms("Signal_ANN_noPT",model_noPT, signal_x_data, signal_pT, signal_CSV, coarse_bins)
	#ANN_Make_Binned_ROC_histograms("BG_ANN_noPT",model_noPT, bg_x_data, bg_pT, bg_CSV, coarse_bins)
	#Make_Binned_ROC_Curves('binned_ANN_noPT','Signal_ANN_noPT','BG_ANN_noPT',coarse_bins, log=True, ANN=True)

	#ANN_Make_Binned_ROC_histograms("Signal_ANN_withPT",model_withPT, signal_x_data, signal_pT, signal_CSV, coarse_bins, addFeature="pT")
	#ANN_Make_Binned_ROC_histograms("BG_ANN_withPT",model_withPT, bg_x_data, bg_pT, bg_CSV, coarse_bins, addFeature="pT")
	#Make_Binned_ROC_Curves('binned_ANN_withPT','Signal_ANN_withPT','BG_ANN_withPT',coarse_bins, log=True, ANN=True)


	#Signal_compare, Background_compare = ANN_x_pT_CSV_label_to_clusterdata(x_data, pT, CSV, labels)

	#efficient_Make_Binned_ROC_histograms('Signal_compare', Signal_compare, cut_bins)    
        #efficient_Make_Binned_ROC_histograms('BG_compare', Background_compare, cut_bins)
	
	#FCM.find_cuts('Thesis_Plots/root_files/BG_compare_histograms.root',cut_bins)
	
	ratio_cuts = [1.833, 1.833, 2.0, 2.0, 2.0, 2.0, 2.167, 2.0, 2.167]

        #efficient_binned_tagged_jets_hist([(Signal_compare, "Signal_compare",(0,2500))],"L4_L1", ratio_cuts, CSV_cuts, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
        #efficient_binned_tagged_jets_hist([(Background_compare, "BG_compare",(0,2500))],"L4_L1", ratio_cuts, CSV_cuts, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	
	#binned_exclusive_tagged_jets_hist("Signal_compare",Signal_compare, "L4/L1", ratio_cuts, CSV_cuts, cut_bins, (0,2500), 60, Difference=False, mode="pT_jet", y_max = 870, Save=True, Stacked=True, AllJets=binned_AllJets_thist)
		
        compare_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_compareL4_L1.root")
        compare_bg_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_compareL4_L1.root")
        binned_ratio_thist = 		FCM.RebinHist(compare_tagged_ratio_file.Get('Signal_compare_Discriminant'),"L4/L1",plot_bins)	#compare_tagged_ratio_file.Get('Signal_compare_Discriminant')		
        binned_bg_ratio_thist = 	FCM.RebinHist(compare_bg_tagged_ratio_file.Get('BG_compare_Discriminant'),"L4/L1",plot_bins)   #compare_bg_tagged_ratio_file.Get('BG_compare_Discriminant')	        

	Efficiency_vs_pT("Signal_ANN_noPT_vs_withPT_vs_L4_L1",[(binned_ANN_noPT_thist,"ANN without p_{T}",3),(binned_ANN_withPT_thist,"ANN with p_{T}",6),(binned_ratio_thist,"L4/L1",2),(binned_CSV_thist,"CSV",4)], binned_AllJets_thist, 0.7, Save=True,legend_shift=True)
        Efficiency_vs_pT("Background_ANN_noPT_vs_withPT_vs_L4_L1",[(binned_bg_ANN_noPT_thist,"ANN without p_{T}",3),(binned_bg_ANN_withPT_thist,"ANN with p_{T}",6),(binned_bg_ratio_thist,"L4/L1",2),(binned_bg_CSV_thist,"CSV",4)], binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True)

	#DrawHistograms([(binned_AllJets_thist,"all jets"), (binned_ANN_noPT_thist, "ANN without p_{T}"),(binned_ANN_withPT_thist, "ANN with p_{T}"),(binned_ratio_thist, "L4/L1"),(binned_CSV_thist,"CSV")], (0,2500), "ANN_noPT_vs_withPT_vs_L4_L1_tagged_jets_vs_pT_binned", 'jet p_{T} (GeV)', "# jets", Save=True,Normalize=False)
		
	
	#Relative_Gain_Plots_ANN()
	'''
	
	#PU study of ANNs

	#print "noPT:"
	#ANN_Make_Binned_ROC_histograms("Signal_ANN",model_noPT, signal_x_data, signal_pT, signal_CSV, cut_bins)
	#ANN_Make_Binned_ROC_histograms("BG_ANN",model_noPT, bg_x_data, bg_pT, bg_CSV, cut_bins)
	#FCM.find_cuts('Thesis_Plots/root_files/BG_ANN_histograms.root',cut_bins,ANN=True)
	#print "withPT:"
	#ANN_Make_Binned_ROC_histograms("Signal_ANN",model_withPT, signal_x_data, signal_pT, signal_CSV, cut_bins,addFeature="pT")
	#ANN_Make_Binned_ROC_histograms("BG_ANN",model_withPT, bg_x_data, bg_pT, bg_CSV, cut_bins,addFeature="pT")
	#FCM.find_cuts('Thesis_Plots/root_files/BG_ANN_histograms.root',cut_bins,ANN=True)
	#print "withPV:"
	#ANN_Make_Binned_ROC_histograms("Signal_ANN",model_withPV, signal_x_data, signal_pT, signal_CSV, cut_bins,addFeature="PV")
	#ANN_Make_Binned_ROC_histograms("BG_ANN",model_withPV, bg_x_data, bg_pT, bg_CSV, cut_bins,addFeature="PV")
	#FCM.find_cuts('Thesis_Plots/root_files/BG_ANN_histograms.root',cut_bins,ANN=True)

	Signal_compare, Background_compare = ANN_x_pT_CSV_label_to_clusterdata(x_data, pT, CSV, labels)
	#efficient_Make_Binned_ROC_histograms('Signal_PU_compare', Signal_compare, cut_bins)    
        #efficient_Make_Binned_ROC_histograms('BG_PU_compare', Background_compare, cut_bins)
	#FCM.find_cuts('Thesis_Plots/root_files/BG_PU_compare_histograms.root',cut_bins)
	
	#ANN_cuts = [0.583, 0.6, 0.617, 0.633, 0.633, 0.667, 0.633, 0.617, 0.717] from model before
	ANN_noPT_cuts = [0.533, 0.55, 0.567, 0.6, 0.617, 0.633, 0.583, 0.617, 0.7]
	ANN_withPT_cuts = [0.6, 0.583, 0.583, 0.583, 0.6, 0.617, 0.617, 0.633, 0.683]
	ANN_withPV_cuts = [0.533, 0.55, 0.583, 0.583, 0.583, 0.683, 0.583, 0.667, 0.75]
	CSV_cuts = [0.65, 0.667, 0.7, 0.7, 0.683, 0.7, 0.667, 0.683, 0.6]
	ratio_cuts = [1.833, 2.0, 2.0, 2.167, 2.0, 2.167, 2.167, 2.167, 3.167]

	#datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_noPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_noPT", (0,2500))]
	#ANN_binned_tagged_jets_hist(datalist, model_noPT, ANN_noPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True)
	#datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_withPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_withPT", (0,2500))]
	#ANN_binned_tagged_jets_hist(datalist, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True,addFeature="pT")
	#datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_withPV", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_withPV", (0,2500))]
	#ANN_binned_tagged_jets_hist(datalist, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True,addFeature="PV")
	
	tagged_ANN_noPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PU_noPTANN.root")
	bg_tagged_ANN_noPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PU_noPTANN.root")
	binned_AllJets_thist = 		FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_AllJets'),"AllJets",plot_bins)   				#tagged_ANN_noPT_file.Get('Signal_PU_noPT_AllJets')		
        binned_CSV_thist = 		FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_CSV'),"CSV",plot_bins)           				#tagged_ANN_noPT_file.Get('Signal_PU_noPT_CSV')		        
        binned_ANN_noPT_thist = 	FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_Discriminant'),"ANN",plot_bins)					#tagged_ANN_noPT_file.Get('Signal_PU_noPT_Discriminant')		
	binned_bg_AllJets_thist = 	FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_AllJets'),"AllJets",plot_bins)    				#bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_AllJets')	        
        binned_bg_CSV_thist = 		FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_CSV'),"CSV",plot_bins)            				#bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_CSV')		        
        binned_bg_ANN_noPT_thist =	FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_Discriminant'),"ANN",plot_bins) 				        #bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_Discriminant')	        
	tagged_ANN_withPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PU_withPTANN.root")
	bg_tagged_ANN_withPT_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PU_withPTANN.root")
	binned_ANN_withPT_thist = 	FCM.RebinHist(tagged_ANN_withPT_file.Get('Signal_PU_withPT_Discriminant'),"ANN",plot_bins)					#tagged_ANN_withPT_file.Get('Signal_PU_withPT_Discriminant')	
	binned_bg_ANN_withPT_thist =	FCM.RebinHist(bg_tagged_ANN_withPT_file.Get('BG_PU_withPT_Discriminant'),"ANN",plot_bins) 					#bg_tagged_ANN_withPT_file.Get('BG_PU_withPT_Discriminant')	
	tagged_ANN_withPV_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PU_withPVANN.root")
	bg_tagged_ANN_withPV_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PU_withPVANN.root")
	binned_ANN_withPV_thist = 	FCM.RebinHist(tagged_ANN_withPV_file.Get('Signal_PU_withPV_Discriminant'),"ANN",plot_bins)					#tagged_ANN_withPV_file.Get('Signal_PU_withPV_Discriminant')	
	binned_bg_ANN_withPV_thist =	FCM.RebinHist(bg_tagged_ANN_withPV_file.Get('BG_PU_withPV_Discriminant'),"ANN",plot_bins) 				        #bg_tagged_ANN_withPV_file.Get('BG_PU_withPV_Discriminant')	

	#efficient_binned_tagged_jets_hist([(Signal_compare, "Signal_PU_compare",(0,2500)),(Background_compare, "BG_PU_compare",(0,2500))],"L4_L1", ratio_cuts, CSV_cuts, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	compare_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_Signal_PU_compareL4_L1.root")
        compare_bg_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BG_PU_compareL4_L1.root")
        binned_ratio_thist = 		FCM.RebinHist(compare_tagged_ratio_file.Get('Signal_PU_compare_Discriminant'),"L4/L1",plot_bins) 	#compare_tagged_ratio_file.Get('Signal_PU_compare_Discriminant')	 
        binned_bg_ratio_thist = 	FCM.RebinHist(compare_bg_tagged_ratio_file.Get('BG_PU_compare_Discriminant'),"L4/L1",plot_bins)        #compare_bg_tagged_ratio_file.Get('BG_PU_compare_Discriminant')		

	#DrawHistograms([(binned_AllJets_thist,"all jets"), (binned_ANN_noPT_thist, "ANN without additional variable"),(binned_ANN_withPT_thist, "ANN with p_{T}"),(binned_ANN_withPV_thist, "ANN with PV"),(binned_ratio_thist, "L4/L1"),(binned_CSV_thist,"CSV")], (0,2500), "ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1_PU_tagged_jets_vs_pT_binned", 'jet p_{T} (GeV)', "# jets", Save=True,Normalize=False)	
	
	#Efficiency_vs_pT("Signal_PU_ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1",[(binned_ANN_noPT_thist,"ANN without p_{T}/PV",3),(binned_ANN_withPT_thist,"ANN with p_{T}",6),(binned_ANN_withPV_thist,"ANN with PV",7),(binned_ratio_thist,"L4/L1",2),(binned_CSV_thist,"CSV",4)], binned_AllJets_thist, 0.7, Save=True,legend_shift=True,LargeLegend=True)
        #Efficiency_vs_pT("Background_PU_ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1",[(binned_bg_ANN_noPT_thist,"ANN without p_{T}/PV",3),(binned_bg_ANN_withPT_thist,"ANN with p_{T}",6),(binned_bg_ANN_withPV_thist,"ANN with PV",7),(binned_bg_ratio_thist,"L4/L1",2),(binned_bg_CSV_thist,"CSV",4)], binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True, LargeLegend=True)


	ANN_efficiency_vs_PU_pT_PV("ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1", signal_x_data, signal_pT, signal_CSV, model_noPT, model_withPT, model_withPV, ANN_noPT_cuts, ANN_withPT_cuts, ANN_withPV_cuts, ratio_cuts, CSV_cuts, cut_bins, 0.7, pT_Cut=200, BG=False)
	ANN_efficiency_vs_PU_pT_PV("ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1_BG", bg_x_data, bg_pT, bg_CSV, model_noPT, model_withPT, model_withPV, ANN_noPT_cuts, ANN_withPT_cuts, ANN_withPV_cuts, ratio_cuts, CSV_cuts, cut_bins, 0.3, pT_Cut=200, BG=True)
	#ANN_efficiency_vs_PU_pT_PV("ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1", signal_x_data, signal_pT, signal_CSV, model_noPT, model_withPT, model_withPV, ANN_noPT_cuts, ANN_withPT_cuts, ANN_withPV_cuts, ratio_cuts, CSV_cuts, cut_bins, 0.7, pT_Cut=1200, BG=False)
	#ANN_efficiency_vs_PU_pT_PV("ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1_BG", bg_x_data, bg_pT, bg_CSV, model_noPT, model_withPT, model_withPV, ANN_noPT_cuts, ANN_withPT_cuts, ANN_withPV_cuts, ratio_cuts, CSV_cuts, cut_bins, 0.3, pT_Cut=1200, BG=True)
	

	
	
	#for appendix: correlation plots between decayvx and pT. also jet-pT vs hadron-pT
		
		
	#also for appendix: differece between hadron and quark

	#QuarkHadronComparison(Signal_2TeV_noPU_String_list, "2TeV", 0.1, 350, BG=False, EarlyBreak=0)
        #QuarkHadronComparison(Signal_4TeV_noPU_String_list, "4TeV", 0.1, 350, BG=False, EarlyBreak=0)
	'''
	quark_file_2TeV = rt.TFile.Open('Thesis_Plots/root_files/QuarkvsHadron_2TeV.root')
	pT_2TeV_hist = quark_file_2TeV.Get("pt_diff"); eta_2TeV_hist = quark_file_2TeV.Get("eta_diff"); phi_2TeV_hist = quark_file_2TeV.Get("phi_diff")
	quark_file_4TeV = rt.TFile.Open('Thesis_Plots/root_files/QuarkvsHadron_4TeV.root')
	pT_4TeV_hist = quark_file_4TeV.Get("pt_diff"); eta_4TeV_hist = quark_file_4TeV.Get("eta_diff"); phi_4TeV_hist = quark_file_4TeV.Get("phi_diff")	
	
	DrawHistograms([(pT_2TeV_hist,"2TeV",3), (pT_4TeV_hist,"4TeV",4)], (0,1), "Quark_vs_Hadron_pT", '#Deltap_{T}/p_{T quark}', "(a.u.)", Save=True,Normalize=True,DrawTitle=False)
	DrawHistograms([(eta_2TeV_hist,"2TeV",3), (eta_4TeV_hist,"4TeV",4)], (0,1), "Quark_vs_Hadron_eta", '#Delta#eta/#eta_{quark}', "(a.u.)", Save=True,Normalize=True,DrawTitle=False)
	DrawHistograms([(phi_2TeV_hist,"2TeV",3), (phi_4TeV_hist,"4TeV",4)], (0,1), "Quark_vs_Hadron_phi", '#Delta#phi/#phi_{quark}', "(a.u.)", Save=True,Normalize=True,DrawTitle=False)
	'''
	'''
	FCM.dR_Dist('2TeV',Signal_2TeV_noPU_String_list, 350, BG=False, EarlyBreak=0)
	FCM.dR_Dist('4TeV',Signal_4TeV_noPU_String_list, 350, BG=False, EarlyBreak=0)
	FCM.dR_Dist('BG',BG_noPU_String_list, 350, BG=False, EarlyBreak=0)
	'''
	'''
	L1 = rt.TH1D("dR_L1","dR_L1,",40,0,0.2)
	L2 = rt.TH1D("dR_L2","dR_L2,",40,0,0.2)
	L3 = rt.TH1D("dR_L3","dR_L3,",40,0,0.2)
	L4 = rt.TH1D("dR_L4","dR_L4,",40,0,0.2)

	dR_dist_Files = []
	for i in range(1,48):
		try:
			dR_dist_Files.append(rt.TFile.Open("Thesis_Plots/root_files/dR_dist_BG_{}.root".format(i)))
			L1.Add(dR_dist_Files[i-1].Get("dR_L1"))
			L2.Add(dR_dist_Files[i-1].Get("dR_L2"))
			L3.Add(dR_dist_Files[i-1].Get("dR_L3"))
			L4.Add(dR_dist_Files[i-1].Get("dR_L4"))
		except:
			print "error: skipped file nr",i
			continue		
	combined_file = rt.TFile("Thesis_Plots/root_files/dR_BG_combined.root","recreate")
	L1.Write()
	L2.Write()
	L3.Write()
	L4.Write()

	print "saved file as Thesis_Plots/root_files/dR_BG_combined.root"
	'''
	'''
	dR_file_2TeV = rt.TFile.Open("Thesis_Plots/root_files/dR_dist_2TeV.root")
	L1_2TeV = dR_file_2TeV.Get("dR_L1"); L2_2TeV = dR_file_2TeV.Get("dR_L2"); L3_2TeV = dR_file_2TeV.Get("dR_L3"); L4_2TeV = dR_file_2TeV.Get("dR_L4")
	dR_file_4TeV = rt.TFile.Open("Thesis_Plots/root_files/dR_dist_4TeV.root")
	L1_4TeV = dR_file_4TeV.Get("dR_L1"); L2_4TeV = dR_file_4TeV.Get("dR_L2"); L3_4TeV = dR_file_4TeV.Get("dR_L3"); L4_4TeV = dR_file_4TeV.Get("dR_L4")
	dR_file_BG = rt.TFile.Open("Thesis_Plots/root_files/dR_BG_combined.root")
	L1_BG = dR_file_BG.Get("dR_L1"); L2_BG = dR_file_BG.Get("dR_L2"); L3_BG = dR_file_BG.Get("dR_L3"); L4_BG = dR_file_BG.Get("dR_L4")
	Draw_dR_dist_Histograms("2TeV", L1_2TeV,L2_2TeV,L3_2TeV,L4_2TeV)
	Draw_dR_dist_Histograms("4TeV", L1_4TeV,L2_4TeV,L3_4TeV,L4_4TeV)
	Draw_dR_dist_Histograms("BG", L1_BG,L2_BG,L3_BG,L4_BG)
	'''







