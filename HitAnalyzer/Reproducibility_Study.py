import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt

def Convert_HA_to_ClusterMatcher(file_path, output_path):
	print "Opening file: "+file_path
	f1 = rt.TFile.Open(file_path)
	tree = f1.Get("demo/tree")
	N = tree.GetEntries()

	NC = np.ndarray((0,26))
        for i in xrange(N):
                if i%1000==0: print "scanning event nr.{}".format(i)
                tree.GetEntry(i)
                for j in xrange(tree.nJets):
                        _NC = np.array([i, tree.jet_bTag[j], 0, tree.jet_pt[j], 0, tree.nClusters_L1004[j],tree.nClusters_L1006[j],tree.nClusters_L1008[j],tree.nClusters_L1010[j],tree.nClusters_L1016[j],tree.nClusters_L2004[j],tree.nClusters_L2006[j],tree.nClusters_L2008[j],tree.nClusters_L2010[j],tree.nClusters_L2016[j],tree.nClusters_L3004[j],tree.nClusters_L3006[j],tree.nClusters_L3008[j],tree.nClusters_L3010[j],tree.nClusters_L3016[j],tree.nClusters_L4004[j],tree.nClusters_L4006[j],tree.nClusters_L4008[j],tree.nClusters_L4010[j],tree.nClusters_L4016[j], tree.nPV])
			NC = np.vstack((NC, _NC))
	
	f1.Close()
	np.save(output_path, NC) 
	print "Saved as: "+output_path

def load_data(file_list):
                data = np.ndarray((0,26))
                for path in file_list:
			try:
                        	temp = np.load(path)
			except IOError:
				temp = np.ndarray((0,26))
                        data = np.vstack((data,temp))
                print 'files loaded'
                return data

def New_Data(title, Signal_data, Background_data):
        n = Signal_data.shape[0]
        if Background_data.shape[0] >1.5*n:
                Background_data = Background_data[:int(1.5*n),:]
        Signal_data[:,0]=1
        Background_data[:,0]=0
	np.random.shuffle(Signal_data)
	np.random.shuffle(Background_data)
        Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
        np.random.shuffle(Data_sample)
        m = Data_sample.shape[0]
        r = int(0.2*m)
	#r = int(m)
        np.save("ANN_data/test_x_{}.npy".format(title),Data_sample[:r,3:])
        np.save("ANN_data/test_y_{}.npy".format(title),Data_sample[:r,0].astype(int))
        np.save("ANN_data/test_CSV_{}.npy".format(title),Data_sample[:r,1])
        np.save("ANN_data/test_pT_{}.npy".format(title),Data_sample[:r,2])
        np.save("ANN_data/train_x_{}.npy".format(title),Data_sample[r:,3:])
        np.save("ANN_data/train_y_{}.npy".format(title),Data_sample[r:,0].astype(int))
        np.save("ANN_data/train_pT_{}.npy".format(title),Data_sample[r:,2])
	np.save("ANN_data/train_CSV_{}.npy".format(title),Data_sample[r:,1])

def Relative_Gain_Plots_ANN():
	from Plots_for_Thesis import Relative_Gains

        ANN_noPT_file = rt.TFile.Open("Thesis_Plots/root_files/repro_Signal_ANN_noPT_tagged_jets_vs_pT_jet_exclusive.root")
        ANN_withPT_file = rt.TFile.Open("Thesis_Plots/root_files/repro_Signal_ANN_withPT_tagged_jets_vs_pT_jet_exclusive.root")
        L4_L1_file = rt.TFile.Open("Thesis_Plots/root_files/repro_Signal_compare_tagged_jets_vs_pT_jet_exclusive_ratio.root")

        noPT_Discriminant_not_CSV = ANN_noPT_file.Get("repro_Signal_ANN_noPT_Discriminant_and_not_CSV")
        noPT_CSV_not_Discriminant = ANN_noPT_file.Get("repro_Signal_ANN_noPT_CSV_not_Discriminant")
        noPT_Discrimant_and_CSV = ANN_noPT_file.Get("repro_Signal_ANN_noPT_Discriminant_and_CSV")
        withPT_Discriminant_not_CSV = ANN_withPT_file.Get("repro_Signal_ANN_withPT_Discriminant_and_not_CSV")
        withPT_CSV_not_Discriminant = ANN_withPT_file.Get("repro_Signal_ANN_withPT_CSV_not_Discriminant")
        withPT_Discrimant_and_CSV = ANN_withPT_file.Get("repro_Signal_ANN_withPT_Discriminant_and_CSV")
        L4_L1_Discriminant_not_CSV = L4_L1_file.Get("repro_Signal_compare_Discriminant_and_not_CSV")
        L4_L1_CSV_not_Discriminant = L4_L1_file.Get("repro_Signal_compare_CSV_not_Discriminant")
        L4_L1_Discrimant_and_CSV = L4_L1_file.Get("repro_Signal_compare_Discriminant_and_CSV")

        ANN_noPT_gains, ANN_noPT_errors, thresholds = Relative_Gains(noPT_Discriminant_not_CSV, noPT_CSV_not_Discriminant, noPT_Discrimant_and_CSV)
        ANN_withPT_gains, ANN_withPT_errors, thresholds = Relative_Gains(withPT_Discriminant_not_CSV, withPT_CSV_not_Discriminant, withPT_Discrimant_and_CSV)
        L4_L1_gains, L4_L1_errors, thresholds = Relative_Gains(L4_L1_Discriminant_not_CSV, L4_L1_CSV_not_Discriminant, L4_L1_Discrimant_and_CSV)

        plt.figure()
        plt.errorbar(thresholds, ANN_noPT_gains, yerr=ANN_noPT_errors, fmt='g', label=r"ANN without $p_T$")
        plt.errorbar(thresholds, ANN_withPT_gains, yerr=ANN_withPT_errors, fmt='magenta', label=r"ANN with $p_T$")
        plt.errorbar(thresholds, L4_L1_gains, yerr=L4_L1_errors, fmt='r', label=r"L4/L1")
        plt.xlim(200,2000)
        plt.ylim(0,2.0)
        plt.xlabel(r"jet $p_T$ threshold (GeV)")
        plt.ylabel(r"gain")
        plt.legend(loc=2)
        plt.savefig("Thesis_Plots/repro_Relative_Gain_Plots.png")
        print "saved figure as Thesis_Plots/repro_Relative_Gain_Plots.png"
        #plt.show()

	import csv
        with open("Thesis_Plots/Relative_Gains_{}.csv".format("repro"), "w") as csvFile:
                writer = csv.writer(csvFile)
                writer.writerow(thresholds)
                writer.writerow(ANN_withPT_gains)
                writer.writerow(ANN_withPT_errors)
        csvFile.close()
        print "saved gains as Thesis_Plots/Relative_Gains_{}.csv".format("repro")

def  pT_PU_from_oldmatcher(data, output_file):
	
        tfile = rt.TFile(output_file+".root","recreate")

        hist_pT = rt.TH1D("pT", "pt", 225, 0, 4500)
        for entry in data[:,3]:
                hist_pT.Fill(entry)
        hist_pT.GetXaxis().SetTitle("jet p_{T} (GeV)")
        hist_pT.GetYaxis().SetTitle('# jets')
        hist_pT.GetYaxis().SetTitleOffset(1.5)
	hist_PU = rt.TH1D("PU", "PU", 100, 0, 100)
	for entry in data[:,25]:
		hist_PU.Fill(entry)
	hist_PU.GetXaxis().SetTitle("# PV")
        hist_PU.GetYaxis().SetTitle('# events')
        hist_PU.GetYaxis().SetTitleOffset(1.5)

        hist_pT.Write()
	hist_PU.Write()
	tfile.Close()

def Plot_MANtag(x_data, pt, model, output_file):
	
	from ANN4Grid import discriminants

	x_Li_Lj=discriminants(x_data)
	x_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
	reduced_pt = pt/200
	
	MANtag = model.predict([x_Li, x_Li_Lj, reduced_pt])

        f2 = rt.TFile(output_file, "RECREATE")
        hist_MANtag = rt.TH1D("MANtag", "MANtag", 100, 0, 1)
        for _MANtag in MANtag:
                hist_MANtag.Fill(_MANtag)
        hist_MANtag.GetXaxis().SetTitle("MANtag")
        hist_MANtag.GetYaxis().SetTitle("# jets")
        hist_MANtag.GetYaxis().SetTitleOffset(1.5)

        hist_MANtag.Write()
        Create_Correlation_Hist("MANtag_vs_pT", "pT", "MANtag", pt, MANtag, (0,3000), (0,1), 100, 100)
        f2.Close()
        print "saved in:", output_file

def Create_Correlation_Hist(title, titleX, titleY, varX, varY, ranX, ranY, binsX, binsY, DrawProfileOnly=False):
        Hist = rt.TH2D(title, title, binsX, ranX[0], ranX[1], binsY, ranY[0], ranY[1])
        rt.gStyle.SetOptStat(0)
        for n,Y  in enumerate(varY):
                Hist.Fill(varX[n], varY[n])
        Hist.GetXaxis().SetTitle(titleX)
        Hist.GetYaxis().SetTitle(titleY)
        Hist.GetYaxis().SetTitleOffset(1.5)
        canvA = rt.TCanvas(title+"_2D", title+"_2D", 600,600)
        rt.gStyle.SetOptStat(0)
        canvA.cd()
        Hist.Draw("COLZ")
	if not DrawProfileOnly: canvA.Write()
        canvB = rt.TCanvas(title+"_Profile", title+"_Profile", 600,600)
        rt.gStyle.SetOptStat(0)
        canvB.cd()
        Hist.ProfileX().Draw("COLZ")
        canvB.Write()
	if not DrawProfileOnly: Hist.ProfileX().Write()

def Inclusive_ROC(output_file, x_data, pt, CSV, labels, model):
	from ANN4Grid import discriminants
	from sklearn.metrics import roc_curve
        x_Li_Lj=discriminants(x_data)
        x_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
        MANtag = model.predict([x_Li, x_Li_Lj, pt/200.])

	fpr, tpr, thresholds = roc_curve(labels, MANtag)
        fpr_csv, tpr_csv, thresholds_csv = roc_curve(labels, CSV)
	
	plt.figure()
        plt.semilogy(tpr,fpr,'r-',label='ANN')
        plt.semilogy(tpr_csv,fpr_csv,'b-',label='CSV')
        plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.ylim(1e-3, 1.)
	plt.xlabel("signal efficiency")
	plt.ylabel("mistag rate")
        plt.legend(loc=4)
	plt.savefig(output_file)
	print "saved as",output_file	

def PU_binned_ROC(output_file, x_data, pt, CSV, labels, model, PU_bins):
	from ANN4Grid import discriminants
	from sklearn.metrics import roc_curve
	if len(PU_bins)<=6:
                colors = ['red','green','blue','orange','brown']
        elif len(PU_bins) <=13:
                colors = ['deepskyblue','rosybrown','olivedrab','royalblue','firebrick','chartreuse','navy','red','darkorchid','lightseagreen','mediumvioletred','blue']
	else:
		print "too many bins, not enough colors!"
		return None

	PU = x_data[:,-1]
        x_Li_Lj=discriminants(x_data)
        x_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
        MANtag = model.predict([x_Li, x_Li_Lj, pt/200.])

	plt.figure()
	for n in range(len(PU_bins)):
		if n == len(PU_bins)-1: continue
		condition =  (PU >= PU_bins[n]) & (PU < PU_bins[n+1])
		fpr, tpr, thresholds = roc_curve(labels[condition], MANtag[condition])
        	fpr_csv, tpr_csv, thresholds_csv = roc_curve(labels[condition], CSV[condition])
	
        	plt.semilogy(tpr,fpr, color=colors[n], linestyle='-', label='MANtag {}_{}'.format(PU_bins[n], PU_bins[n+1]))
        	plt.semilogy(tpr_csv,fpr_csv, color=colors[n], linestyle='--', label='CSV {}_{}'.format(PU_bins[n], PU_bins[n+1]))

        plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.ylim(1e-3, 1.)
	plt.xlabel("signal efficiency")
	plt.ylabel("mistag rate")
        plt.legend(loc=4)
	plt.savefig(output_file)
	print "saved as",output_file	

def Efficiency_vs_PU(output_file, x_data, pT, CSV, model, ANN_Cuts, CSV_Cuts, bins, y_max, pT_Cut=200, BG=False):
	assert x_data.shape[1]==21, "x_data does not contain PV. Make sure it is made from a PU sample and has shape (x, 21)."
	assert x_data.shape[0] == len(pT) == len(CSV), "data inputs need to have the same length"

	from Plots_for_Thesis import ANN_bin_selection, ANN_functional_shape
	ran = (0,80)
        nbins = 80
        import array
        if BG:
                bins_ = array.array('d',[0.0, 11.0]+range(19,41,8)+[42.0,  52.0, 80])
        else:
                bins_ = array.array('d',[0.0, 11.0]+range(15,41,4)+[42.0, 52.0, 58.0, 65.0, 80])

        if pT_Cut >= 1200:
                bins_ = array.array('d',[0.0, 20.0, 40.0, 80.0])

	AllJets_Hist = rt.TH1D("AllJets","AllJets",nbins,ran[0],ran[1])
        ANN_Hist = rt.TH1D("MANtag","MANtag",nbins,ran[0],ran[1])
        CSV_Hist = rt.TH1D("CSV","CSV",nbins,ran[0],ran[1])
        AllJets_Hist = AllJets_Hist.Rebin(len(bins_)-1,"AllJets",bins_)
        ANN_Hist = ANN_Hist.Rebin(len(bins_)-1,"MANtag",bins_)
        CSV_Hist = CSV_Hist.Rebin(len(bins_)-1,"CSV",bins_)

        pred_y = model.predict(ANN_functional_shape(x_data)+[pT/200.])
        bin_numbers = ANN_bin_selection(pT,bins)

	for i,pT_value in enumerate(pT):
                        if pT_value < pT_Cut: continue
                        if bin_numbers[i] == -100: continue
                        AllJets_Hist.Fill(x_data[i,-1])
                        if CSV[i] >= CSV_Cuts[bin_numbers[i]]: CSV_Hist.Fill(x_data[i,-1])
                        if pred_y[i] >= ANN_Cuts[bin_numbers[i]]: ANN_Hist.Fill(x_data[i,-1])

	canvas = rt.TCanvas('canvas','canvas',600,600)
        rt.gStyle.SetOptTitle(0)
        legend = rt.TLegend(0.1,0.9,0.35,0.75)
        ANN_Graph = rt.TGraphAsymmErrors()
        CSV_Graph = rt.TGraphAsymmErrors()
        ANN_Graph.Divide(ANN_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        CSV_Graph.Divide(CSV_Hist,AllJets_Hist,"cl=0.683 b(1,1) mode")
        ANN_Graph.SetLineColor(6)
        CSV_Graph.SetLineColor(4)
        legend.AddEntry(ANN_Graph, "MANtag", "LEP")
        legend.AddEntry(CSV_Graph, "CSV", "LEP")
        if BG:
                ANN_Graph.GetYaxis().SetTitle('mistag rate')
        else:
                ANN_Graph.GetYaxis().SetTitle('efficiency')
        ANN_Graph.GetYaxis().SetTitleOffset(1.5)
        ANN_Graph.SetMinimum(0.)
        ANN_Graph.SetMaximum(y_max)
        ANN_Graph.Draw()
        CSV_Graph.Draw("SAME")
        legend.Draw()
        canvas.SaveAs(output_file)

def Input_Reconstruction(x_data, pt, model, output_file, MANtag_range=(0,1)):
	from ANN4Grid import discriminants
	assert x_data.shape[1]==21, "x_data has wrong shape"
	pt.shape = (pt.size, 1)

	#predict from model
	x_Li_Lj=discriminants(x_data)
        x_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
        MANtag = model.predict([x_Li, x_Li_Lj, pt/200.])

	#define MANtag constraint
	condition =  (MANtag >= MANtag_range[0]) & (MANtag <= MANtag_range[1])

	#separate out input variables
	Li_list = []
	for i in range(0,20):
		Li_list.append(x_data[:,i])
		Li_list[i].shape = (Li_list[i].size, 1)
	PU = x_data[:,20] #not really input variable but still interesting to check correlation
	PU.shape = (PU.size,1)

	#plot output
	hist_MANtag = rt.TH1D("MANtag_selected", "MANtag_selected", 100, 0, 1)
	for _MANtag in MANtag[condition]:
                hist_MANtag.Fill(_MANtag)
	hist_MANtag.GetXaxis().SetTitle("MANtag")
        hist_MANtag.GetYaxis().SetTitle("# jets")
        hist_MANtag.GetYaxis().SetTitleOffset(1.5)
	hist_MANtag.SetLineColor(4)
	hist_MANtag_full = rt.TH1D("MANtag_full", "MANtag_full", 100, 0, 1)
	for _MANtag in MANtag:
                hist_MANtag_full.Fill(_MANtag)
	hist_MANtag_full.GetXaxis().SetTitle("MANtag")
        hist_MANtag_full.GetYaxis().SetTitle("# jets")
        hist_MANtag_full.GetYaxis().SetTitleOffset(1.5)
	hist_MANtag_full.SetLineColor(2)

	#draw output on canvas
	canv_MANtag = rt.TCanvas("MANtag", "MANtag", 600, 600)
	rt.gStyle.SetOptStat(0)	
	legend = rt.TLegend(0.65,0.75, 0.9, 0.9)
	legend.AddEntry(hist_MANtag, "{} #leq MANtag #leq {}".format(MANtag_range[0],  MANtag_range[1]), "LEP")
	legend.AddEntry(hist_MANtag_full, "full", "LEP")
	canv_MANtag.cd()
	hist_MANtag_full.Draw()
	hist_MANtag.Draw("SAME")
	legend.Draw()

	#save to file and add input variables
	f2 = rt.TFile(output_file, "RECREATE")
	canv_MANtag.Write()
	#hist_MANtag.Write()
	Add_MANtag_Input_Histogram("pT", (0,3000), 100, pt, MANtag, condition, legend)
	Add_MANtag_Input_Histogram("PU", (0,80), 80, PU, MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L1_004", (0,30), 30, Li_list[0], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L1_006", (0,35), 35, Li_list[1], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L1_008", (0,41), 41, Li_list[2], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L1_010", (0,42), 42, Li_list[3], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L1_016", (0,51), 51, Li_list[4], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L2_004", (0,51), 51, Li_list[5], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L2_006", (0,66), 66, Li_list[6], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L2_008", (0,84), 84, Li_list[7], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L2_010", (0,92), 92, Li_list[8], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L2_016", (0,91), 91, Li_list[9], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L3_004", (0,42), 42, Li_list[10], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L3_006", (0,48), 48, Li_list[11], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L3_008", (0,46), 46, Li_list[12], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L3_010", (0,48), 48, Li_list[13], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L3_016", (0,58), 58, Li_list[14], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L4_004", (0,42), 42, Li_list[15], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L4_006", (0,49), 49, Li_list[16], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L4_008", (0,53), 53, Li_list[17], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L4_010", (0,54), 54, Li_list[18], MANtag, condition, legend)
	Add_MANtag_Input_Histogram("L4_016", (0,61), 61, Li_list[19], MANtag, condition, legend)
	f2.Close()
	print "saved in:", output_file

def Add_MANtag_Input_Histogram(var_title, var_range, nbins, variable, MANtag, condition, legend):
	hist = rt.TH1D(var_title+"_selected", var_title+"_selected", nbins, var_range[0], var_range[1])
        for _var in variable[condition]:
                hist.Fill(_var)
        hist.GetXaxis().SetTitle(var_title)
        hist.GetYaxis().SetTitle("(a.u.)")
        hist.GetYaxis().SetTitleOffset(1.5)
	hist.SetLineColor(4)
	hist_full = rt.TH1D(var_title+"_full", var_title+"_full", nbins, var_range[0], var_range[1])
        for _var in variable:
                hist_full.Fill(_var)
        hist_full.GetXaxis().SetTitle(var_title)
        hist_full.GetYaxis().SetTitle("(a.u.)")
        hist_full.GetYaxis().SetTitleOffset(1.5)
	hist_full.SetLineColor(2)

	canv = rt.TCanvas(var_title, var_title, 600, 600)
        rt.gStyle.SetOptStat(0)
	rt.gStyle.SetOptTitle(0)		
        canv.cd()
        hist.DrawNormalized()
        hist_full.DrawNormalized("SAME")
        legend.Draw()
	canv.Write()
	#hist.Write()
	Create_Correlation_Hist("MANtag_vs_"+var_title, var_title, "MANtag", variable, MANtag, var_range, (0,1), nbins, 100, DrawProfileOnly=True)


if __name__ == "__main__":

	signal_path = "/eos/user/m/msommerh/HitAnalyzer/{}/M{}/flatTuple_{}.root"
	bg_path = "/eos/user/m/msommerh/HitAnalyzer/QCD/{}/flatTuple_{}.root"

	output_directory = "/afs/cern.ch/work/m/msommerh/public/ClusterMatcher"

        M0_list = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']
        bin_list = [('170','300'), ('300','470'), ('470','600'), ('600','800'), ('800','1000'), ('1000','1400'), ('1400','1800'), ('1800','2400'), ('2400','3200'), ('3200', 'Inf')]


	##Make ClusterMatcher files out of matched HitAnalyzer files
        #for M0 in M0_list[5:10]:
        #        for i in range(1,16):
        #                Convert_HA_to_ClusterMatcher(signal_path.format('2017',M0,i), output_directory+"/2017/M{}/CM_{}.npy".format(M0, i))
        #                Convert_HA_to_ClusterMatcher(signal_path.format('2018',M0,i), output_directory+"/2018/M{}/CM_{}.npy".format(M0, i))
        #for bin_ in bin_list:
        #        for i in range(1,42):
        #                Convert_HA_to_ClusterMatcher(bg_path.format(bin_[0]+"to"+bin_[1],i), output_directory+"/QCD/{}/CM_{}.npy".format(bin_[0]+"to"+bin_[1], i))


	##Create Datasets for ANNs from new HitAnalyzer files
	#signal_list = []
	#bg_list = []
	#signal_CM_path = output_directory+"/{}/M{}/CM_{}.npy"
	#bg_CM_path = output_directory+"/QCD/{}/CM_{}.npy"
        #for M0 in M0_list[5:10]:
        #        for i in range(1,16):
	#		signal_list.append(signal_CM_path.format(2017, M0, i))
       	#		signal_list.append(signal_CM_path.format(2018, M0, i))
	#signal_data = load_data(signal_list)
	#for bin_ in bin_list:
        #        for i in range(1,42):
       	#		bg_list.append(bg_CM_path.format(bin_[0]+"to"+bin_[1], i))
	#bg_data = load_data(bg_list)
	#New_Data("repro", signal_data, bg_data)

	#Create Datasets for ANNs from old HitAnalyzer files
	#signal_list = []
	#bg_list = []
	#CM_path = "/afs/cern.ch/work/m/msommerh/public/repro_matched_clusters/{}/MatchedClusters_{}_{}_{}.npy"
        #for M0 in M0_list[5:]:
        #        for i in range(1,31):
	#		signal_list.append(CM_path.format("M"+M0, 2017, "M"+M0, i))
       	#		signal_list.append(CM_path.format("M"+M0, 2018, "M"+M0, i))
	#signal_data = load_data(signal_list)
	#for bin_ in bin_list:
        #        for i in range(1,85):
       	#		bg_list.append(CM_path.format(bin_[0]+"to"+bin_[1], "QCD", bin_[0]+"to"+bin_[1], i))
	#bg_data = load_data(bg_list)
	##New_Data("repro", signal_data, bg_data)

	#
        ##generate weight function for all cases:
	#from FinalClusterMatcher import find_weight_function
        #print 'weitghts with pt'
        #find_weight_function(bg_data, signal_data, 'jet_pT', 5)
        #print 'weights with PV'
        #find_weight_function(bg_data, signal_data, 'PU',1)


	##Plot pT and PU new data matched with old framework
	#signal_list = []
        #bg_list = []
        #CM_path = "/afs/cern.ch/work/m/msommerh/public/repro_matched_clusters/{}/MatchedClusters_{}_{}_{}.npy"
        #for M0 in M0_list[5:]:
        #        for i in range(1,31):
        #               signal_list.append(CM_path.format("M"+M0, 2017, "M"+M0, i))
        #               signal_list.append(CM_path.format("M"+M0, 2018, "M"+M0, i))
        #signal_data = load_data(signal_list)
        #for bin_ in bin_list:
        #        for i in range(1,85):
        #               bg_list.append(CM_path.format(bin_[0]+"to"+bin_[1], "QCD", bin_[0]+"to"+bin_[1], i))
        #bg_data = load_data(bg_list)
	#pT_PU_from_oldmatcher(signal_data, "pT_PU_newdata_oldmatcher_signal")
	#pT_PU_from_oldmatcher(bg_data, "pT_PU_newdata_oldmatcher_bg")

	#import FinalClusterMatcher as FCM
	#Signal_4TeV_PU = FCM.load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
        #Signal_2TeV_PU = FCM.load_data('2TeV-Signal_PU','matched_clusters/Signal_PU/',19)
        #Signal_both_PU = np.vstack((Signal_4TeV_PU,Signal_2TeV_PU))
        #Background_PU = FCM.load_data('BG_PU','matched_clusters/BG_PU/',499)
 	#pT_PU_from_oldmatcher(Signal_both_PU, "pT_PU_olddata_oldmatcher_signal")
	#pT_PU_from_oldmatcher(Background_PU, "pT_PU_olddata_oldmatcher_bg")
       
	
	import keras as kr
        model_noPT = kr.models.load_model("Submitted_Models/repro/model_withPU_functional.h5")
        model_withPT = kr.models.load_model("Submitted_Models/repro/model_withPU_functional_withPT.h5")
        #model_withPV = kr.models.load_model("Submitted_Models/model_withPU_functional_withPV.h5")
        model_type = "functional"
        print "model loaded"

        x_data = np.load("ANN_data/repro2/test_x.npy")
        CSV = np.load("ANN_data/repro2/test_CSV.npy")
        labels = np.load("ANN_data/repro2/test_y.npy")
        pT = np.load("ANN_data/repro2/test_feature.npy")
        print "data loaded"
	#x_data = np.load("ANN_data/test_x_repro.npy")
        #CSV = np.load("ANN_data/test_CSV_repro.npy")
        #labels = np.load("ANN_data/test_y_repro.npy")
        #pT = np.load("ANN_data/test_pT_repro.npy")
        #print "data loaded"

        signal_x_data = x_data[labels==1]
        signal_pT = pT[labels==1]
        signal_CSV = CSV[labels==1]
        bg_x_data = x_data[labels==0]
        bg_pT = pT[labels==0]
        bg_CSV = CSV[labels==0]	
	

	#Plot_MANtag(signal_x_data, signal_pT, model_withPT, "repro_plots/new_old_signal_MANtag_dist.root")
	#Plot_MANtag(bg_x_data, bg_pT, model_withPT, "repro_plots/new_old_bg_MANtag_dist.root")
	'''	
	import FinalClusterMatcher as FCM
	from repro_Plots_for_Thesis import ANN_Make_Binned_ROC_histograms, ANN_x_pT_CSV_label_to_clusterdata, efficient_Make_Binned_ROC_histograms, Efficiency_vs_pT, ANN_binned_tagged_jets_hist, efficient_binned_tagged_jets_hist, find_cuts
	
	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500] #used for the pT-depending cuts
	plot_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500] #used for a nice looking plot
	
	##noPT:
        ANN_Make_Binned_ROC_histograms("Signal_ANN_repro",model_noPT, signal_x_data, signal_pT, signal_CSV, cut_bins)
        ANN_Make_Binned_ROC_histograms("BG_ANN_repro",model_noPT, bg_x_data, bg_pT, bg_CSV, cut_bins)
        ANN_noPT_cuts, CSV_cuts = find_cuts('Thesis_Plots/root_files/BG_ANN_repro_histograms.root',cut_bins,ANN=True)
        ##withPT:
        ANN_Make_Binned_ROC_histograms("Signal_ANN_repro",model_withPT, signal_x_data, signal_pT, signal_CSV, cut_bins,addFeature="pT")
        ANN_Make_Binned_ROC_histograms("BG_ANN_repro",model_withPT, bg_x_data, bg_pT, bg_CSV, cut_bins,addFeature="pT")
        ANN_withPT_cuts, CSV_cuts = find_cuts('Thesis_Plots/root_files/BG_ANN_repro_histograms.root',cut_bins,ANN=True)
        ##withPV:
        #ANN_Make_Binned_ROC_histograms("Signal_ANN_repro",model_withPV, signal_x_data, signal_pT, signal_CSV, cut_bins,addFeature="PV")
        #ANN_Make_Binned_ROC_histograms("BG_ANN_repro",model_withPV, bg_x_data, bg_pT, bg_CSV, cut_bins,addFeature="PV")
        #ANN_withPV_cuts, CSV_cuts = find_cuts('Thesis_Plots/root_files/BG_ANN_repro_histograms.root',cut_bins,ANN=True)

	
        Signal_compare, Background_compare = ANN_x_pT_CSV_label_to_clusterdata(x_data, pT, CSV, labels)
        efficient_Make_Binned_ROC_histograms('Signal_PU_compare_repro', Signal_compare, cut_bins)
        efficient_Make_Binned_ROC_histograms('BG_PU_compare_repro', Background_compare, cut_bins)
        delta_cuts, ratio_cuts, CSV_cuts = find_cuts('Thesis_Plots/root_files/BG_PU_compare_repro_histograms.root',cut_bins)

        #ANN_noPT_cuts = [0.866666666667, 0.833333333333, 0.833333333333, 0.816666666667, 0.783333333333, 0.7, 0.716666666667]
        #ANN_withPT_cuts = [0.916666666667, 0.9, 0.883333333333, 0.866666666667, 0.816666666667, 0.7, 0.45]
        #ANN_withPV_cuts = [0.7, 0.7, 0.716666666667, 0.716666666667, 0.733333333333, 0.7, 0.283333333333]
        #CSV_cuts = [0.7, 0.683333333333, 0.7, 0.7, 0.716666666667, 0.7, 0.2833333333333]
        #ratio_cuts = [1.0, 1.0, 1.33333333333, 1.33333333333, 1.33333333333, 1.83333333333, 0.166666666667]
	
        datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_noPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_noPT", (0,2500))]
        ANN_binned_tagged_jets_hist(datalist, model_noPT, ANN_noPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True)
        datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_withPT", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_withPT", (0,2500))]
        ANN_binned_tagged_jets_hist(datalist, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True,addFeature="pT")
        #datalist = [(signal_x_data, signal_pT, signal_CSV, "Signal_PU_withPV", (0,2500)),(bg_x_data, bg_pT, bg_CSV, "BG_PU_withPV", (0,2500))]
        #ANN_binned_tagged_jets_hist(datalist, model_withPT, ANN_withPV_cuts, CSV_cuts, cut_bins, 60, mode="pT_jet",Save=True,addFeature="PV")
	
        tagged_ANN_noPT_file =  rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_Signal_PU_noPTANN.root")
        bg_tagged_ANN_noPT_file =       rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_BG_PU_noPTANN.root")
        binned_AllJets_thist =          FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_AllJets'),"AllJets",plot_bins)         #tagged_ANN_noPT_file.Get('Signal_PU_noPT_AllJets')     	                          
        binned_CSV_thist =              FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_CSV'),"CSV",plot_bins)                 #tagged_ANN_noPT_file.Get('Signal_PU_noPT_CSV')         	                          
        binned_ANN_noPT_thist =         FCM.RebinHist(tagged_ANN_noPT_file.Get('Signal_PU_noPT_Discriminant'),"ANN",plot_bins)        #tagged_ANN_noPT_file.Get('Signal_PU_noPT_Discriminant')	                          
        binned_bg_AllJets_thist =       FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_AllJets'),"AllJets",plot_bins)          #bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_AllJets')      	                          
        binned_bg_CSV_thist =           FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_CSV'),"CSV",plot_bins)                  #bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_CSV')          	                          
        binned_bg_ANN_noPT_thist =      FCM.RebinHist(bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_Discriminant'),"ANN",plot_bins)         #bg_tagged_ANN_noPT_file.Get('BG_PU_noPT_Discriminant') 	                          
        tagged_ANN_withPT_file =        rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_Signal_PU_withPTANN.root")
        bg_tagged_ANN_withPT_file =     rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_BG_PU_withPTANN.root")
        binned_ANN_withPT_thist =       FCM.RebinHist(tagged_ANN_withPT_file.Get('Signal_PU_withPT_Discriminant'),"ANN",plot_bins)      #tagged_ANN_withPT_file.Get('Signal_PU_withPT_Discriminant')	                                
        binned_bg_ANN_withPT_thist =    FCM.RebinHist(bg_tagged_ANN_withPT_file.Get('BG_PU_withPT_Discriminant'),"ANN",plot_bins)       #bg_tagged_ANN_withPT_file.Get('BG_PU_withPT_Discriminant') 	                                
        #tagged_ANN_withPV_file =        rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_Signal_PU_withPVANN.root")
        #bg_tagged_ANN_withPV_file =     rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_BG_PU_withPVANN.root")
        #binned_ANN_withPV_thist =       FCM.RebinHist(tagged_ANN_withPV_file.Get('Signal_PU_withPV_Discriminant'),"ANN",plot_bins)        #tagged_ANN_withPV_file.Get('Signal_PU_withPV_Discriminant')	                              
        #binned_bg_ANN_withPV_thist =    FCM.RebinHist(bg_tagged_ANN_withPV_file.Get('BG_PU_withPV_Discriminant'),"ANN",plot_bins)         #bg_tagged_ANN_withPV_file.Get('BG_PU_withPV_Discriminant') 	                              

        efficient_binned_tagged_jets_hist([(Signal_compare, "Signal_PU_compare",(0,2500)),(Background_compare, "BG_PU_compare",(0,2500))],"L4_L1", ratio_cuts, CSV_cuts, cut_bins, 60, Difference=False, mode="pT_jet",Save=True)
        compare_tagged_ratio_file =     rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_Signal_PU_compareL4_L1.root")
        compare_bg_tagged_ratio_file =  rt.TFile.Open("Thesis_Plots/root_files/repro_tagged_jets_vs_pT_jet_BG_PU_compareL4_L1.root")
        binned_ratio_thist =            FCM.RebinHist(compare_tagged_ratio_file.Get('Signal_PU_compare_Discriminant'),"L4/L1",plot_bins)      #compare_tagged_ratio_file.Get('Signal_PU_compare_Discriminant')	 
        binned_bg_ratio_thist =         FCM.RebinHist(compare_bg_tagged_ratio_file.Get('BG_PU_compare_Discriminant'),"L4/L1",plot_bins)       #compare_bg_tagged_ratio_file.Get('BG_PU_compare_Discriminant') 	 

        #DrawHistograms([(binned_AllJets_thist,"all jets"), (binned_ANN_noPT_thist, "ANN without additional variable"),(binned_ANN_withPT_thist, "ANN with p_{T}"),(binned_ANN_withPV_thist, "ANN with PV"),(binned_ratio_thist, "L4/L1"),(binned_CSV_thist,"CSV")], (0,2500), "ANN_noPT_vs_withPT_vs_withPV_vs_L4_L1_PU_tagged_jets_vs_pT_binned", 'jet p_{T} (GeV)', "# jets", Save=True,Normalize=False)

	Efficiency_vs_pT("repro_Signal_PU_ANN_noPT_vs_withPT_vs_L4_L1",[(binned_ANN_noPT_thist,"ANN without p_{T}",3),(binned_ANN_withPT_thist,"ANN with p_{T}",6),(binned_ratio_thist,"L4/L1",2),(binned_CSV_thist,"CSV",4)], binned_AllJets_thist, 0.7, Save=True,legend_shift=False,LargeLegend=False)
        Efficiency_vs_pT("repro_Background_PU_ANN_noPT_vs_withPT_vs_L4_L1",[(binned_bg_ANN_noPT_thist,"ANN without p_{T}",3),(binned_bg_ANN_withPT_thist,"ANN with p_{T}",6),(binned_bg_ratio_thist,"L4/L1",2),(binned_bg_CSV_thist,"CSV",4)], binned_bg_AllJets_thist, 0.3, Save=True,legend_shift=False, BG=True, LargeLegend=True)


	from Plots_for_Thesis import ANN_exclusive_tagged_jets_hist, binned_exclusive_tagged_jets_hist

        binned_AllJets_thist = tagged_ANN_noPT_file.Get('Signal_PU_noPT_AllJets')     	                          

     	ANN_exclusive_tagged_jets_hist('repro_Signal_ANN_noPT', model_noPT, signal_x_data, signal_pT, signal_CSV, ANN_noPT_cuts, CSV_cuts, cut_bins, (0,2500), 60, mode="pT_jet", y_max=5000, Save=True, Stacked=True, AllJets=binned_AllJets_thist )
        ANN_exclusive_tagged_jets_hist('repro_Signal_ANN_withPT', model_withPT, signal_x_data, signal_pT, signal_CSV, ANN_withPT_cuts, CSV_cuts, cut_bins, (0,2500), 60, mode="pT_jet", y_max=5000, Save=True,addFeature="pT", Stacked=True, AllJets=binned_AllJets_thist ) 

	binned_exclusive_tagged_jets_hist("repro_Signal_compare",Signal_compare, "L4/L1", ratio_cuts, CSV_cuts, cut_bins, (0,2500), 60, Difference=False, mode="pT_jet", y_max = 5000, Save=True, Stacked=True, AllJets=binned_AllJets_thist)

	Relative_Gain_Plots_ANN()
	'''	
	
	
	##new_old ROC curves:
	from Plots_for_Thesis import ANN_Make_Binned_ROC_histograms, Make_Binned_ROC_Curves
	coarse_bins = [0,1200,1800,2500]
	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500]
	PU_bins = [0,15,30,45, 90]
	ANN_withPT_cuts = [0.65000000000000002, 0.6333333333333333, 0.6333333333333333, 0.6166666666666667, 0.6333333333333333, 0.6333333333333333, 0.6166666666666667, 0.6166666666666667]
	CSV_cuts = [0.66666666666666663, 0.66666666666666663, 0.71666666666666667, 0.71666666666666667, 0.71666666666666667, 0.73333333333333328, 0.75, 0.78333333333333333]
	#ANN_Make_Binned_ROC_histograms("repro_new_old_Signal_ANN",model_withPT, signal_x_data, signal_pT, signal_CSV, coarse_bins, addFeature="pT")
        #ANN_Make_Binned_ROC_histograms("repro_new_old_BG_ANN",model_withPT, bg_x_data, bg_pT, bg_CSV, coarse_bins, addFeature="pT")
        #Make_Binned_ROC_Curves('repro_new_old_ANN','repro_new_old_Signal_ANN','repro_new_old_BG_ANN',coarse_bins, log=True, ANN=True)
	#Inclusive_ROC("repro_plots/inclusive_ROC_new_old.png", x_data, pT, CSV, labels, model_withPT)
	#PU_binned_ROC("repro_plots/PU_binned_ROC_new_old.png", x_data, pT, CSV, labels, model_withPT, PU_bins)
	#Efficiency_vs_PU("repro_plots/Eff_vs_PU_new_old.png", signal_x_data, signal_pT, signal_CSV, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 0.7, pT_Cut=200, BG=False)
	#Efficiency_vs_PU("repro_plots/Mistag_vs_PU_new_old.png", bg_x_data, bg_pT, bg_CSV, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 0.3, pT_Cut=200, BG=True)
	
	Input_Reconstruction(signal_x_data, signal_pT, model_withPT, "repro_plots/input_reconstruction_new_old_signal.root", MANtag_range=(0.62,0.68))
	Input_Reconstruction(bg_x_data, bg_pT, model_withPT, "repro_plots/input_reconstruction_new_old_bg.root", MANtag_range=(0.62,0.68))

	'''
	##old_old ROC curves:
	import keras as kr
        model_noPT = kr.models.load_model("Submitted_Models/model_withPU_functional.h5")
        model_withPT = kr.models.load_model("Submitted_Models/model_withPU_functional_withPT.h5")
        #model_withPV = kr.models.load_model("Submitted_Models/model_withPU_functional_withPV.h5")
        model_type = "functional"
        print "model loaded"

	x_data = np.load("Submitted_Models/data/withPU_both_withPT/test_x.npy")
        CSV = np.load("Submitted_Models/data/withPU_both_withPT/test_CSV.npy")
        labels = np.load("Submitted_Models/data/withPU_both_withPT/test_y.npy")
        pT = np.load("Submitted_Models/data/withPU_both_withPT/test_feature.npy")
        print "data loaded"

        signal_x_data = x_data[labels==1]
        signal_pT = pT[labels==1]
        signal_CSV = CSV[labels==1]
        bg_x_data = x_data[labels==0]
        bg_pT = pT[labels==0]
        bg_CSV = CSV[labels==0]

	#Plot_MANtag(signal_x_data, signal_pT, model_withPT, "repro_plots/old_old_signal_MANtag_dist.root")
        #Plot_MANtag(bg_x_data, bg_pT, model_withPT, "repro_plots/old_old_bg_MANtag_dist.root")

	from Plots_for_Thesis import ANN_Make_Binned_ROC_histograms, Make_Binned_ROC_Curves
	coarse_bins = [0,1200,1800,2500]
	cut_bins = [0, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500]
	ANN_withPT_cuts = [0.6, 0.583, 0.583, 0.583, 0.6, 0.617, 0.617, 0.633, 0.683]
        CSV_cuts = [0.65, 0.667, 0.7, 0.7, 0.683, 0.7, 0.667, 0.683, 0.6]
	#ANN_Make_Binned_ROC_histograms("repro_old_old_Signal_ANN",model_withPT, signal_x_data, signal_pT, signal_CSV, coarse_bins, addFeature="pT")
        #ANN_Make_Binned_ROC_histograms("repro_old_old_BG_ANN",model_withPT, bg_x_data, bg_pT, bg_CSV, coarse_bins, addFeature="pT")
        #Make_Binned_ROC_Curves('repro_old_old_ANN','repro_old_old_Signal_ANN','repro_old_old_BG_ANN',coarse_bins, log=True, ANN=True)
	#Inclusive_ROC("repro_plots/inclusive_ROC_old_old.png", x_data, pT, CSV, labels, model_withPT)
	#PU_binned_ROC("repro_plots/PU_binned_ROC_old_old.png", x_data, pT, CSV, labels, model_withPT, PU_bins)
	#Efficiency_vs_PU("repro_plots/Eff_vs_PU_old_old.png", signal_x_data, signal_pT, signal_CSV, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 0.7, pT_Cut=200, BG=False)
	#Efficiency_vs_PU("repro_plots/Mistag_vs_PU_old_old.png", bg_x_data, bg_pT, bg_CSV, model_withPT, ANN_withPT_cuts, CSV_cuts, cut_bins, 0.3, pT_Cut=200, BG=True)
	'''
	

