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
                if particle[4] >= minR and particle[3] >= minPT2:
                        L1_th += particle[dR_tag]
                        L2_th += particle[dR_tag+1]
                        L3_th += particle[dR_tag+2]
                        L4_th += particle[dR_tag+3]

        fig2, ax2 = plt.subplots(1,2,figsize=(9,5))
        #fig2.suptitle('Hit Clusters per Layer inside dR<'+str(dR)+' on '+title+' sample')
        ax2[0].bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
	ax2[0].set_title("jet-pT>{}GeV/c".format(minPT1))
        ax2[0].set_ylabel('Clusters')
        ax2[0].set_xticks([0.5,1.5,2.5,3.5])
        ax2[0].set_xticklabels(['L1','L2','L3','L4'])
        ax2[1].bar([0.5,1.5,2.5,3.5],[L1_th, L2_th, L3_th, L4_th],align='center')
	ax2[1].set_title("jet-pT>{}GeV/c".format(minPT2))
        #ax2[1].set_ylabel('[a.u.]')
        ax2[1].set_xticks([0.5,1.5,2.5,3.5])
        ax2[1].set_xticklabels(['L1','L2','L3','L4'])
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

def General_Make_ROC_Curves(title, histlist,log=False, print_cut=False):
	'''histlist entry: (signal_hist, bg_hist, title, diff( True or False), signal_ZeroDiv, bg_ZeroDiv,linestyle(optional))'''
        #hsv = plt.get_cmap('hsv')
        #color = hsv(np.linspace(0,1.0,len(bins)-1))
        #color = ['b', 'g', 'r', 'c', 'm', 'y']
        if len(histlist)<=6:
                color = ['red','green','blue','orange','brown','black']
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
                plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"$\epsilon$_background")
                plt.legend(loc=4)
        else:
                #plt.plot([0,0],[0,0],'k-',label = 'L4/L1')
                #plt.plot([0,0],[0,0],'k-.',label = 'CSV')
                #plt.plot([0,1],[0.9,0.9],'k:',label="10% mistag")
                plt.plot([0,1],[0.9,0.9],'k:')
                plt.xlabel(r"$\epsilon$_signal")
                plt.ylabel(r"1-$\epsilon$_background")
                plt.legend(loc=3)
        plt.title(title+"_ROC-Curves")

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
	plt.plot(Cuts,signal_eff,'.',label=r'$\epsilon$ signal')
	plt.plot(Cuts,1-bg_eff,'.',label=r'1-$\epsilon$ background')
        plt.xlabel(r"$\epsilon$ signal")
        plt.ylabel(r"1-$\epsilon$ background")
        plt.legend(loc=5)
        plt.savefig("Thesis_Plots/efficiency_vs_cut_{}.png".format(title))
        print "saved as Thesis_Plots/efficiency_vs_cut_{}.png".format(title)

def Efficiency_vs_pT(title,histlist, hist_all_jets,y_max,Save=False,legend_shift=False,BG=False):
        """plots for each histogram of tagged jets given in a list of tuples (histogram, title) the efficiency for each bin, where the x-axis corresponds to the feature given as string (see FeatureDict)."""
        canvas = rt.TCanvas('canvas','canvas',600,600)
        if legend_shift:
                legend = rt.TLegend(0.1,0.1,0.35,0.25)
        else:
                legend = rt.TLegend(0.1,0.9,0.35,0.75)
        graphlist = []
        for n,hist in enumerate(histlist):
                graphlist.append(rt.TGraphAsymmErrors())
                #if n==0: graphlist[n].SetTitle(title+"_vs_jet-pT")
                graphlist[n].Divide(hist[0],hist_all_jets,"cl=0.683 b(1,1) mode")
                legend.AddEntry(graphlist[n], histlist[n][1],"LEP")
                graphlist[n].SetLineColor(n+2)
                if n<1:
                        graphlist[n].GetXaxis().SetTitle("jet-pT (GeV/c)")
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

def Make_Binned_ROC_Curves(title,Signal_title,Background_title,bins, diff=False,log=False):
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

        Signal_ZeroDiv = np.loadtxt("Thesis_Plots/root_files/{}_ZeroDiv.csv".format(Signal_title),delimiter=',')
        Signal_file = rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(Signal_title),"READ")
        Background_ZeroDiv = np.loadtxt("Thesis_Plots/root_files/{}_ZeroDiv.csv".format(Background_title),delimiter=',')
        Background_file =    rt.TFile("Thesis_Plots/root_files/{}_histograms.root".format(Background_title),"READ")

        plt.figure("ROC")
        plt.clf()

        for bin_ in range(len(bins)-1):
                Dis_Signal_Eff = FCM.Get_ROC_Efficiencies(Signal_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),ran,nbins,Signal_ZeroDiv[bin_])
                Dis_BG_Eff = FCM.Get_ROC_Efficiencies(Background_file.Get(dis_string+str(bins[bin_])+"_"+str(bins[bin_+1])),ran,nbins,Background_ZeroDiv[bin_])
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
                AllJetsHistlist[n].GetXaxis().SetTitle(mode)
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
                try:
                        L_R = particle[16]/float(particle[13])
                        if L_R >= Ratio_Cuts[bin_number]: Ratio_Hist.Fill(particle[25])
                except ZeroDivisionError:
                        continue

        AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        Delta_Hist = Delta_Hist.Rebin(len(bins_)-1,"Delta",bins_)
        Ratio_Hist = Ratio_Hist.Rebin(len(bins_)-1,"Ratio",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)

        #Make Graphs and draw them
        canvas = rt.TCanvas('canvas','canvas',600,600)
        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        Delta_Graph = rt.TGraphAsymmErrors()
        Ratio_Graph = rt.TGraphAsymmErrors()
        CSV_Graph = rt.TGraphAsymmErrors()
        Ratio_Graph.SetTitle(title+"_vs_PU_pT{}{}".format(pT_Mode,pT_Cut))
        Delta_Graph.Divide(Delta_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Ratio_Graph.Divide(Ratio_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        Delta_Graph.SetLineColor(3)
        Ratio_Graph.SetLineColor(2)
        CSV_Graph.SetLineColor(4)
        legend.AddEntry(Delta_Graph, "L4-L1", "LEP")
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
        Delta_Graph.Draw("SAME")
        CSV_Graph.Draw("SAME")
        legend.Draw()
        canvas.SaveAs('Thesis_Plots/'+title+"_vs_PU_pT{}{}.png".format(pT_Mode,pT_Cut))

def DrawHistograms(Histograms, ran, title, xlabel, ylabel, Save=False,Normalize=True, t_sleep=0):
        """Draws multiple histograms neatly on a canvas when given to the function as list of tuples (root histogram, title)."""
        canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
	histlist = []
        if len(Histograms) > 1:
                rt.gStyle.SetOptStat(0)#something is wrong with this
                legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for nr, Histogram in enumerate(Histograms):
		histlist.append(Histogram[0])
		if nr < 3:
                	histlist[nr].SetLineColor(nr+2)
		else:
			histlist[nr].SetLineColor(nr+3)
                if nr == 0:
			histlist[nr].SetTitle(title)
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
	
	'''
	Signal_4TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_4TeV-Signal.npy')
        Signal_2TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_2TeV-Signal.npy')
        Signal_both_noPU = np.vstack((Signal_4TeV_noPU,Signal_2TeV_noPU))
        Background_noPU = FCM.load_data('BG','matched_clusters/BG_noPU/',31, old_version=True)#31
        '''
		
	Signal_4TeV_PU = FCM.load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
        Signal_2TeV_PU = FCM.load_data('2TeV-Signal_PU','matched_clusters/Signal_PU/',19)
        Signal_both_PU = np.vstack((Signal_4TeV_PU,Signal_2TeV_PU))
        Background_PU = FCM.load_data('BG_PU','matched_clusters/BG_PU/',499)	
	
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
	Histograms = [(jet_pT_hist_BG,"Background"),(jet_pT_hist_2TeV,"2TeV-Signal"),(jet_pT_hist_4TeV,"4TeV-Signal")]
	DrawHistograms(Histograms, (0,2500), "jet_pT", "jet_pT (GeV/c)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''
	'''
	Efficient_Sample_Analysis("Signal_2TeV_decayvx", Signal_2TeV_noPU, (0,35), "decayvx")
	Efficient_Sample_Analysis("Signal_4TeV_decayvx", Signal_4TeV_noPU, (0,35), "decayvx")
	Efficient_Sample_Analysis("BG_decayvx", Background_noPU, (0,0.2), "decayvx")
	tfile1 = rt.TFile.Open("Thesis_Plots/root_files/Signal_2TeV_decayvx.root")
	decayvx_hist_2TeV = tfile1.Get("decay_vertex_R")
	tfile2 = rt.TFile.Open("Thesis_Plots/root_files/Signal_4TeV_decayvx.root")
	decayvx_hist_4TeV = tfile2.Get("decay_vertex_R")
	tfile3 = rt.TFile.Open("Thesis_Plots/root_files/BG_decayvx.root")
	decayvx_hist_BG = tfile3.Get("decay_vertex_R")
	Histograms = [(decayvx_hist_2TeV,"2TeV-Signal"),(decayvx_hist_4TeV,"4TeV-Signal")]
	DrawHistograms(Histograms, (0,35), "decay_vertex_R_Signal", "decay vertex R (cm)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	DrawHistograms([(decayvx_hist_BG,"Background")], (0,0.2), "decay_vertex_R_Background", "decay vertex R (cm)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''
	'''
	Efficient_Sample_Analysis("Signal_2TeV_CSV", Signal_2TeV_noPU, (0,1), "CSV")
	Efficient_Sample_Analysis("Signal_4TeV_CSV", Signal_4TeV_noPU, (0,1), "CSV")
	Efficient_Sample_Analysis("BG_CSV", Background_noPU, (0,1), "CSV")
	tfile1 = rt.TFile.Open("Thesis_Plots/root_files/Signal_2TeV_CSV.root")
	CSV_hist_2TeV = tfile1.Get("CSV")
	tfile2 = rt.TFile.Open("Thesis_Plots/root_files/Signal_4TeV_CSV.root")
	CSV_hist_4TeV = tfile2.Get("CSV")
	tfile3 = rt.TFile.Open("Thesis_Plots/root_files/BG_CSV.root")
	CSV_hist_BG = tfile3.Get("CSV")
	Histograms = [(CSV_hist_BG,"Background"),(CSV_hist_2TeV,"2TeV-Signal"),(CSV_hist_4TeV,"4TeV-Signal")]
	DrawHistograms(Histograms, (0,1), "CSV_b-tag", "p(b-jet)","(a.u.)", Save=True,Normalize=True, t_sleep=0)
	'''

	#nice plot illustrating cluster matching

	'''	
	with open("Grid.pkl",) as f:   #open coordinates of DetUnits for visual reference
               Grid = pickle.load(f)
	ax = FCM.Initialize3DPlot('Cluster Matching', 'x', 'y', 'z', grid=Grid)
	#HitClusters = CM.ClusterMatch(rt.TFile.Open(Signal_4TeV_noPU_String.format(20)), 0.1, 350, HadronsNotQuarks=True, Plot=True, Axes=ax, Save=False, dR_dist=False, LayerHist=False, EarlyBreak=5)
	'''
	#satisfactory example of decay after first layer falsely reconstructed by CSV	
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
	FCM.Efficient_SeparateLayerHist([(Signal_2TeV_noPU,"2TeV"),(Signal_4TeV_noPU,"4TeV"),(Background_noPU,"Background")], (0,30), 0.1 , minPT=200,jet_pT=True, Save=True)
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

	diff_histlist = [(signal_4TeV_diff_21, bg_diff_21, "L2-L1", "diff", 0, 0), (signal_4TeV_diff_31, bg_diff_31, "L3-L1", "diff", 0, 0), (signal_4TeV_diff_41, bg_diff_41, "L4-L1", "diff", 0, 0), (signal_4TeV_diff_32, bg_diff_32, "L3-L2", "diff", 0, 0), (signal_4TeV_diff_42, bg_diff_42, "L4-L2", "diff", 0, 0), (signal_4TeV_diff_43, bg_diff_43, "L4-L3", "diff", 0, 0)]

	General_Make_ROC_Curves("ratio_taggers_dR{}".format(dR), ratio_histlist,log=False)
	General_Make_ROC_Curves("diff_taggers_dR{}".format(dR), diff_histlist,log=False)
	'''
	#Discriminant Histograms corresponding to the best discriminants: L4/L1 at dR<0.1 and L4-L1 at dR<0.04
	'''
	tfile_2TeV_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_2TeV_dR{}.root".format(0.04))
	tfile_4TeV_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_dR{}.root".format(0.04))
	tfile_BG_004 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_dR{}.root".format(0.04))
	tfile_2TeV_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_2TeV_dR{}.root".format(0.1))
	tfile_4TeV_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_4TeV_dR{}.root".format(0.1))
	tfile_BG_01 = rt.TFile.Open("Thesis_Plots/root_files/GeneralROCHists_BG_dR{}.root".format(0.1))

	signal_2TeV_ratio = tfile_2TeV_01.Get("L4_L1")
	signal_4TeV_ratio = tfile_4TeV_01.Get("L4_L1")
	bg_ratio = tfile_BG_01.Get("L4_L1")
	signal_2TeV_diff = tfile_2TeV_004.Get("L4-L1")
	signal_4TeV_diff = tfile_4TeV_004.Get("L4-L1")
	bg_diff = tfile_BG_004.Get("L4-L1")
	signal_2TeV_CSV = tfile_2TeV_01.Get("CSV")	
	signal_4TeV_CSV = tfile_4TeV_01.Get("CSV")	
	bg_CSV = tfile_BG_01.Get("CSV")	

	#DrawHistograms([(bg_ratio,'Background'),(signal_2TeV_ratio,'2TeV'),(signal_4TeV_ratio,'4TeV')], (0,10), 'L4_L1_Discriminant_Hist', 'L4/L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	#DrawHistograms([(bg_diff,'Background'),(signal_2TeV_diff,'2TeV'),(signal_4TeV_diff,'4TeV')], (-22,22), 'L4-L1_Discriminant_Hist', 'L4-L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	#DrawHistograms([(bg_CSV,'Background'),(signal_2TeV_CSV,'2TeV'),(signal_4TeV_CSV,'4TeV')], (0,1), 'CSV_Discriminant_Hist', 'p(b-jet)', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	'''
	
	#best ROC-curves global and single pT threshold (discriminant variables from above)
	
	'''	
	bg_ZeroDiv_21, bg_ZeroDiv_31, bg_ZeroDiv_41, bg_ZeroDiv_32, bg_ZeroDiv_42, bg_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_BG_dR{}.csv".format(0.1),delimiter=',')
	signal_4TeV_ZeroDiv_21, signal_4TeV_ZeroDiv_31, signal_4TeV_ZeroDiv_41, signal_4TeV_ZeroDiv_32, signal_4TeV_ZeroDiv_42, signal_4TeV_ZeroDiv_43 = np.loadtxt("Thesis_Plots/root_files/GeneralROCHists_ZeroDiv_4TeV_dR{}.csv".format(0.1),delimiter=',')
	'''
	'''
	#best_histlist = [(signal_4TeV_ratio, bg_ratio, "L4/L1", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41), (signal_4TeV_diff, bg_diff, "L4-L1", "diff", 0,0),(signal_4TeV_CSV, bg_CSV, "CSV", "CSV", 0,0)]
	#General_Make_ROC_Curves("Best_Taggers", best_histlist,log=False,print_cut=True)
	
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

	diff_HPT_histlist = [(signal_4TeV_diff, bg_diff, "L4-L1 ($p_T$>200GeV)", "diff", 0,0, 'blue','-'),(signal_4TeV_diff_HPT, bg_diff_HPT, "L4-L1 ($p_T$>1200GeV)", "diff", 0,0, 'blue','--'),(signal_4TeV_CSV, bg_CSV, "CSV ($p_T$>200GeV)", "CSV", 0,0,'red','-'),(signal_4TeV_CSV_HPT, bg_CSV_HPT, "CSV ($p_T$>1200GeV)", "CSV", 0,0,'red','--')]
	ratio_HPT_histlist = [(signal_4TeV_ratio, bg_ratio, "L4/L1 ($p_T$>200GeV)", "ratio", signal_4TeV_ZeroDiv_41, bg_ZeroDiv_41, 'blue','-'),(signal_4TeV_ratio_HPT, bg_ratio_HPT, "L4/L1 ($p_T$>1200GeV)", "ratio", signal_4TeV_ZeroDiv_41_HPT, bg_ZeroDiv_41_HPT, 'blue', '--'),(signal_4TeV_CSV, bg_CSV, "CSV ($p_T$>200GeV)", "CSV", 0,0, 'red','-'),(signal_4TeV_CSV_HPT, bg_CSV_HPT, "CSV ($p_T$>1200GeV)", "CSV", 0,0, 'red', '--')]
	
	#General_Make_ROC_Curves("diff_Taggers_high_pT", diff_HPT_histlist,log=False,print_cut=True)
	#General_Make_ROC_Curves("ratio_Taggers_high_pT", ratio_HPT_histlist,log=False,print_cut=True)
	'''
	
	#Discriminant histograms above pT threshold
	'''
	DrawHistograms([(bg_ratio_HPT,'Background'),(signal_4TeV_ratio_HPT,'4TeV')], (0,10), 'L4_L1_Discriminant_Hist_high_pT', 'L4/L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	DrawHistograms([(bg_diff_HPT,'Background'),(signal_4TeV_diff_HPT,'4TeV')], (-22,22), 'L4-L1_Discriminant_Hist_high_pT', 'L4-L1', '(a.u)', Save=True, Normalize=True, t_sleep=0)
	DrawHistograms([(bg_CSV_HPT,'Background'),(signal_4TeV_CSV_HPT,'4TeV')], (0,1), 'CSV_Discriminant_Hist_high_pT', 'p(b-jet)', '(a.u)', Save=True, Normalize=True, t_sleep=0)
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
	
	#study on pT-bins

	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500] #used for the pT-depending cuts
	coarse_bins = [0,1200,1800,2500]	
	
	#efficient_Make_Binned_ROC_histograms('Signal', Signal_both_noPU, cut_bins)    
        #efficient_Make_Binned_ROC_histograms('BG', Background_noPU, cut_bins)
	
	#FCM.find_cuts('Thesis_Plots/root_files/BG_binned_histograms.root',cut_bins)

	delta_cuts_noPU = [3, 5, 6, 6, 7, 8, 8, 9, 9]
        ratio_cuts_noPU = [1.833, 1.833, 2.0, 2.0, 2.0, 2.167, 2.167, 2.167, 2.167]
        CSV_cuts_noPU = [0.633, 0.65, 0.667, 0.683, 0.7, 0.717, 0.733, 0.767, 0.767]

	#efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,2500))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        #efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,2500))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	#efficient_binned_tagged_jets_hist([(Background_noPU, "BG",(0,2500))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        #efficient_binned_tagged_jets_hist([(Background_noPU, "BG",(0,2500))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	'''
	binned_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_SignalL4-L1.root")
        binned_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_SignalL4_L1.root")
	binned_bg_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BGL4-L1.root")
        binned_bg_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_pT_jet_BGL4_L1.root")
	binned_AllJets_thist = 		binned_tagged_diff_file.Get('Signal_AllJets')		#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_AllJets'),"AllJets",plot_bins)
        binned_CSV_thist = 		binned_tagged_diff_file.Get('Signal_CSV')		#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_CSV'),"CSV",plot_bins)
        binned_diff_thist = 		binned_tagged_diff_file.Get('Signal_Discriminant')	#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_Discriminant'),"L4-L1",plot_bins)
        binned_ratio_thist = 		binned_tagged_ratio_file.Get('Signal_Discriminant')	#FCM.RebinHist(binned_tagged_ratio_file.Get('Signal_Discriminant'),"L4/L1",plot_bins)
	binned_bg_AllJets_thist = 	binned_bg_tagged_diff_file.Get('BG_AllJets')		#FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_AllJets'),"AllJets",plot_bins)
        binned_bg_CSV_thist = 		binned_bg_tagged_diff_file.Get('BG_CSV')		#FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_CSV'),"CSV",plot_bins)
        binned_bg_diff_thist = 		binned_bg_tagged_diff_file.Get('BG_Discriminant')	#FCM.RebinHist(binned_bg_tagged_diff_file.Get('BG_Discriminant'),"L4-L1",plot_bins)
        binned_bg_ratio_thist = 	binned_bg_tagged_ratio_file.Get('BG_Discriminant')	#FCM.RebinHist(binned_bg_tagged_ratio_file.Get('BG_Discriminant'),"L4/L1",plot_bins)
	'''
	'''
	Efficiency_vs_pT("Signal_binned",[(binned_ratio_thist,"L4/L1"),(binned_diff_thist,"L4-L1"),(binned_CSV_thist,"CSV")], binned_AllJets_thist, 0.6, Save=True,legend_shift=True)
        Efficiency_vs_pT("Background_binned",[(binned_bg_ratio_thist,"L4/L1"),(binned_bg_diff_thist,"L4-L1"),(binned_bg_CSV_thist,"CSV")], binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True)
	'''
	#DrawHistograms([(binned_AllJets_thist,"all jets"), (binned_ratio_thist, "L4/L1"), (binned_diff_thist, "L4-L1"), (binned_CSV_thist,"CSV")], (0,2500), "tagged_jets_vs_pT_binned", 'jet-pT', "# jets", Save=True,Normalize=False)
	'''
	efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,30))],"L4-L1", delta_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=True, mode="decay_vx",Save=True)
        efficient_binned_tagged_jets_hist([(Signal_both_noPU, "Signal",(0,30))],"L4_L1", ratio_cuts_noPU, CSV_cuts_noPU, cut_bins, 60, Difference=False, mode="decay_vx",Save=True)
	
	binned_dvx_tagged_diff_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_decay_vx_SignalL4-L1.root")
        binned_dvx_tagged_ratio_file = 	rt.TFile.Open("Thesis_Plots/root_files/binned_tagged_jets_vs_decay_vx_SignalL4_L1.root")
	binned_dvx_AllJets_thist = 		binned_dvx_tagged_diff_file.Get('Signal_AllJets')	#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_AllJets'),"AllJets",plot_bins)
        binned_dvx_CSV_thist = 			binned_dvx_tagged_diff_file.Get('Signal_CSV')		#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_CSV'),"CSV",plot_bins)
        binned_dvx_diff_thist = 		binned_dvx_tagged_diff_file.Get('Signal_Discriminant')	#FCM.RebinHist(binned_tagged_diff_file.Get('Signal_Discriminant'),"L4-L1",plot_bins)
        binned_dvx_ratio_thist = 		binned_dvx_tagged_ratio_file.Get('Signal_Discriminant')	#FCM.RebinHist(binned_tagged_ratio_file.Get('Signal_Discriminant'),"L4/L1",plot_bins)

	DrawHistograms([(binned_dvx_AllJets_thist,"all jets"), (binned_dvx_ratio_thist, "L4/L1"), (binned_dvx_diff_thist, "L4-L1"), (binned_dvx_CSV_thist,"CSV")], (0,30), "tagged_jets_vs_decayvx_binned", 'decay vertex R', "# jets", Save=True,Normalize=False)
	'''
	'''
	#efficient_Make_Binned_ROC_histograms('Signal_coarse-binned', Signal_both_noPU, coarse_bins)
        #efficient_Make_Binned_ROC_histograms('BG_coarse-binned', Background_noPU, coarse_bins)
	Make_Binned_ROC_Curves('coarse_binned_ratio','Signal_coarse-binned','BG_coarse-binned',coarse_bins, diff=False,log=False)
	Make_Binned_ROC_Curves('coarse_binned_diff','Signal_coarse-binned','BG_coarse-binned',coarse_bins, diff=True,log=False)
	'''

	#PU study of cut-based taggers
	'''
	#PU_histograms(Signal_2TeV_PU, Signal_4TeV_PU, Background_PU)	
	
	PU_file = rt.TFile.Open("Thesis_Plots/root_files/PU_distributions.root")
	PU_2TeV = PU_file.Get("PU_2TeV"); PU_4TeV = PU_file.Get("PU_4TeV"); PU_BG = PU_file.Get("PU_BG")
	DrawHistograms([(PU_BG,"Background"), (PU_2TeV,"2TeV Signal"), (PU_4TeV,"4TeV Signal")], (0,80), "PU_distributions", '#PV', "(a.u.)", Save=True,Normalize=True)
	'''
	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500] #used for the pT-depending cuts
	'''
	efficient_Make_Binned_ROC_histograms('Signal_PU_binned', Signal_both_PU, cut_bins)    
        efficient_Make_Binned_ROC_histograms('BG_PU_binned', Background_PU, cut_bins)

	FCM.find_cuts('Thesis_Plots/root_files/BG_PU_binned_histograms.root',cut_bins)
	'''
	delta_cuts_PU = [4, 5, 6, 7, 8, 8, 9, 9, 11]
        ratio_cuts_PU = [1.833, 2.0, 2.0, 2.167, 2.167, 2.167, 2.167, 2.167, 2.167]
        CSV_cuts_PU = [0.65, 0.667, 0.683, 0.7, 0.717, 0.733, 0.75, 0.783, 0.783]	
	'''
	efficient_binned_tagged_jets_hist([(Signal_both_PU, "Signal_PU",(0,2500))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        efficient_binned_tagged_jets_hist([(Signal_both_PU, "Signal_PU",(0,2500))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	efficient_binned_tagged_jets_hist([(Background_PU, "BG_PU",(0,2500))],"L4-L1", delta_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=True, mode="pT_jet",Save=True)
        efficient_binned_tagged_jets_hist([(Background_PU, "BG_PU",(0,2500))],"L4_L1", ratio_cuts_PU, CSV_cuts_PU, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
	
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

	binned_efficiency_vs_PU('Efficiency', Signal_both_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.8, pT_Cut=200, pT_Mode="jet")
        binned_efficiency_vs_PU('Efficiency', Signal_both_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.8, pT_Cut=1200, pT_Mode="jet")

        binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.3, pT_Cut=200, pT_Mode="jet",BG=True)
        binned_efficiency_vs_PU('Mistag-Rate', Background_PU, delta_cuts_PU, ratio_cuts_PU, CSV_cuts_PU, cut_bins, 0.3, pT_Cut=1200, pT_Mode="jet",BG=True)


	#comparison of different ANN models


	#PU study of ANNs


	#for appendix: correlation plots between decayvx and pT. also jet-pT vs hadron-pT
		
		
	#also for appendix: differece between hadron and quark


        
