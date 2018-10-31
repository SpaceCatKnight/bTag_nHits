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

FeatureDict = {"nEvent":0, "nParticle":1, "nLayer":2, "nModule":3, "nCluster":4, "pT_hadron":5, "pdgId":6, "decay_vx":7, "dR":8, "bTag":9, "pT_jet":10}

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
	"""takes the path toa rot file as input and returns the global cartesian coordinates of all the detUnits hit until event denoted by EarlyBreak and saves them to a pkl-file. Used for visual reference in 3D-plots"""
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
	#c1 = rt.TCanvas('canvas','canvas',600,600)
	#rt.gStyle.SetOptStat(0)
        #hist.GetXaxis().SetTitle("R")
        #hist.GetYaxis().SetTitle("# det_units")
        #hist.GetYaxis().SetTitleOffset(1.5)
        #hist.Draw()
	#c1.SaveAs("det_unit_R.png")

def PlotTrajectory(vx,p,ax,T_max,res,col,lwidth,lstyle):
	"""plots the trajectory of a particle starting at vertex point vx=(x,y,z) with velocity/momentum vector p=(px,py,pz) on 3D plot ax with length T_max, resolution res, color col, line width lwidth and line style lstyle"""
        Tx,Ty,Tz=[],[],[]
        v_t = normalize(np.array([p[0],p[1],p[2]]))
        for t in np.linspace(0,T_max,res):
                Tx.append(vx[0]+t*v_t[0])
                Ty.append(vx[1]+t*v_t[1])
                Tz.append(vx[2]+t*v_t[2])
        ax.plot(xs=Tx,ys=Ty,zs=Tz,color=col,linewidth=lwidth,linestyle=lstyle)#plots all the tracks

def DiscriminatorHist(title,discriminator,Bdata,Backgrounddata,bins,ran,xlabel):
	"""creates histogram on signal data Bdata and Backgrounddata according to a given discriminator function"""
	HistB = rt.TH1D(title,title,bins,ran[0],ran[1])
	HistBackground = rt.TH1D(title,title,bins,ran[0],ran[1])
	for particle in Bdata:
		HistB.Fill(discriminator(particle))
	for particle in Backgrounddata:
		HistBackground.Fill(discriminator(particle))
	canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        HistB.GetXaxis().SetTitle(xlabel)
        HistB.GetYaxis().SetTitle('[a.u.]')
        HistB.GetYaxis().SetTitleOffset(1.5)
        HistB.SetLineColor(2)
        HistBackground.SetLineColor(3)
        HistB.DrawNormalized()
        HistBackground.DrawNormalized('SAME')
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        legend.AddEntry(HistB,'B-hadrons')
        legend.AddEntry(HistBackground,'Background')
        legend.Draw()
        canvas.SaveAs('discriminants/'+title+'.png')

def Li_Lj_Hist1D(i, j, Bdatalist, Backgrounddata, bins, ran, dR, Difference=False, Abs=False, dR_check=False, pT_hadron=0,pT_jet=0, Save=False):
	"""creates a histogram of the discriminant Li-Lj (Difference = True) or Li/Lj (Difference = False) on signal data put in a list of tuples (data,title) if dR_ckeck==True, the dR value is used as a threshold on the matched clusters. pT_hadron and pT_jet can be used as additional thresholds."""
        pT_string = ''
	if pT_hadron>0: pT_string = pT_string+'_pTH_'+str(pT_hadron)
	if pT_jet>0: pT_string = pT_string+'_pTJ_'+str(pT_jet)
	if Difference:
                if Abs:
                        title = "abs_L"+str(i)+"-L"+str(j)+"_dR_"+str(dR)+pT_string
                        xlabel = "|L"+str(i)+"-L"+str(j)+"|"
                        ran = (0,25)
                else:
                        title = "L"+str(i)+"-L"+str(j)+"_dR_"+str(dR)+pT_string
                        xlabel = "L"+str(i)+"-L"+str(j)
                        ran = (-25,25)
                bins = ran[1]-ran[0]
        else:
                title = "L"+str(i)+"_L"+str(j)+"_dR_"+str(dR)+pT_string
                xlabel = "L"+str(i)+"/L"+str(j)
        BG_zero_div = 0
        BG_else = 0
        HistBackground = rt.TH1D("Background",title,bins,ran[0],ran[1])
        HistBlist = []
        B_zero_div_list = []
        B_else_list = []
        for n,Bdata in enumerate(Bdatalist):
                print "working on",Bdata[1]
                B_zero_div = 0
                B_else = 0
		ntot = 0.
		npass = 0.
                HistBlist.append(rt.TH1D(Bdata[1],title,bins,ran[0],ran[1]))
                HistBlist[n].SetLineColor(n+3)
                for particle in Bdata[0]:
			ntot +=1.
			if (particle[0][0][5] < pT_hadron) or (particle[0][0][10] < pT_jet): continue
                        npass += 1.
			if Difference:
				if Abs:
                                	L = abs(Li_Lj_diff(i,j,particle, dR, dR_check=dR_check))
                        	else:
                                	L = Li_Lj_diff(i,j,particle, dR, dR_check=dR_check)
                        else:
                                L = Li_Lj_ratio(i,j,particle, dR, dR_check=dR_check)
                        	B_else += 1
                       	if L == 'zero-div':
                               	B_zero_div += 1
                        else:
                                HistBlist[n].Fill(L)
		print npass,"of",ntot,"passed the pT threshold ->",npass/ntot*100,"%"
                B_zero_div_list.append(B_zero_div)
                B_else_list.append(B_else)
        print "working on background"
        for particle in Backgrounddata:
		if (particle[0][0][5] < pT_hadron) or (particle[0][0][10] < pT_jet): continue
                if Difference:
                        if Abs:
                                L = abs(Li_Lj_diff(i,j,particle, dR, dR_check=dR_check))
                        else:
                                L = Li_Lj_diff(i,j,particle, dR, dR_check=dR_check)
                else:
                        L = Li_Lj_ratio(i,j,particle, dR, dR_check=dR_check)
                	BG_else += 1
                if L == 'zero-div':
                       	BG_zero_div += 1
                else:
                        HistBackground.Fill(L)
        canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        HistBlist[0].GetXaxis().SetTitle(xlabel)
        HistBlist[0].GetYaxis().SetTitle('[a.u.]')
        HistBlist[0].GetYaxis().SetTitleOffset(1.5)
        HistBackground.SetLineColor(2)
        HistBlist[0].DrawNormalized()
        if len(HistBlist)>1:
                for n,HistB in enumerate(HistBlist):
                        if n>0: HistB.DrawNormalized('SAME')
        HistBackground.DrawNormalized('SAME')
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for n,HistB in enumerate(HistBlist):
                legend.AddEntry(HistB,Bdatalist[n][1])
        legend.AddEntry(HistBackground,'Background')
        legend.Draw()
        if BG_zero_div != 0 and BG_else != 0: print "Background sample has",BG_zero_div,"out of",BG_else,"zero-divsion entries ->",100*BG_zero_div/float(BG_else),"%"
        for n,Bdata in enumerate(Bdatalist):
                if B_zero_div_list[n] != 0 and B_else_list[n] != 0: print Bdata[1],"sample has",B_zero_div_list[n],"out of",B_else_list[n],"zero-divsion entries ->",100*B_zero_div_list[n]/float(B_else_list[n]),"%"
        if Save:
                Tfile = rt.TFile("histogram_files/"+title+"histograms.root","recreate")
                for Hist in HistBlist: Hist.Write()
                HistBackground.Write()
		if not Difference:
			zero_div = {"background": BG_zero_div}
			for n,Bdata in enumerate(Bdatalist):
				zero_div[Bdata[1]] = B_zero_div_list[n]
                	with open("histogram_files/"+title+"_zero_div.pkl", 'w') as f:
                       		pickle.dump(zero_div, f)
                canvas.SaveAs('discriminants/'+title+'.png')
                print 'saved as discriminants/'+title+'.png'

def CSV_Hist(Bdatalist, Backgrounddata, bins, pT_hadron=0,pT_jet=0, Save=False):
        ran = (0,1)
	pT_string = ''
	if pT_hadron>0: pT_string = pT_string+'_pTH_'+str(pT_hadron)
	if pT_jet>0: pT_string = pT_string+'_pTJ_'+str(pT_jet)
        title = "CSV"+pT_string
        xlabel = "CSV"
        HistBackground = rt.TH1D("Background",title,bins,ran[0],ran[1])
        HistBlist = []
        for n,Bdata in enumerate(Bdatalist):
                print "working on",Bdata[1]
                HistBlist.append(rt.TH1D(Bdata[1],title,bins,ran[0],ran[1]))
                HistBlist[n].SetLineColor(n+3)
                for particle in Bdata[0]:
			if (particle[0][0][5] < pT_hadron) or (particle[0][0][10] < pT_jet): continue
                        HistBlist[n].Fill(particle[0][0][9])
        print "working on background"
        for particle in Backgrounddata:
                HistBackground.Fill(particle[0][0][9])
        canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        HistBlist[0].GetXaxis().SetTitle(xlabel)
        HistBlist[0].GetYaxis().SetTitle('[a.u.]')
        HistBlist[0].GetYaxis().SetTitleOffset(1.5)
        HistBackground.SetLineColor(2)
        HistBlist[0].DrawNormalized()
        if len(HistBlist)>1:
                for n,HistB in enumerate(HistBlist):
                        if n>0: HistB.DrawNormalized('SAME')
        HistBackground.DrawNormalized('SAME')
        legend = rt.TLegend(0.9,0.9,0.65,0.75)
        for n,HistB in enumerate(HistBlist):
                legend.AddEntry(HistB,Bdatalist[n][1])
        legend.AddEntry(HistBackground,'Background')
        legend.Draw()
        if Save:
                Tfile = rt.TFile("histogram_files/"+title+"histograms.root","recreate")
                for Hist in HistBlist: Hist.Write()
                HistBackground.Write()
                canvas.SaveAs('discriminants/'+title+'.png')
                print 'saved as discriminants/'+title+'.png'

def Li_Lj_Hist2D(title,i,j,data,ran,dR,Save=False):
	"""creates a 2D histogram of the amount of pixel clusters in Li vs Lj. dR here is only used for the title. All the clusters are taken into account"""
	bins = ran[1] - ran[0]
	Hist = rt.TH2I(title,title,bins,ran[0],ran[1],bins,ran[0],ran[1])
	for particle in data:
		Li,Lj = 0,0
		for cluster in particle:
			if cluster[0][2] == i: Li += 1
			if cluster[0][2] == j: Lj += 1
		Hist.Fill(Li,Lj)
	canvas = rt.TCanvas('canvas','canvas',600,600)
        canvas.SetTitle(title)
        rt.gStyle.SetOptStat(0)
        Hist.GetXaxis().SetTitle('L'+str(i))
        Hist.GetYaxis().SetTitle('L'+str(j))
        Hist.GetYaxis().SetTitleOffset(1.5)
        Hist.Draw("COLZ")
        if Save: 
		canvas.SaveAs('discriminants/L'+str(i)+'_L'+str(j)+'Hist_2D_DR'+str(dR)+title+'.png')
		print 'saved as discriminants/L'+str(i)+'_L'+str(j)+'Hist_2D_DR'+str(dR)+title+'.png'
	sleep(5)

def Correlation_Hist2D(title, data, feature_x, feature_y, ran_x, ran_y, bins_x, bins_y, Profiling=False, Save=False):
	"""Creates a 2D histogram of features saved in the cluster matching function (see FeatureDict)"""
	f_x, f_y = FeatureDict[feature_x], FeatureDict[feature_y]
	Hist = rt.TH2I(title,title,bins_x,ran_x[0],ran_x[1],bins_y,ran_y[0],ran_y[1])
	for particle in data:
		Hist.Fill(particle[0][0][f_x], particle[0][0][f_y])
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


def Li_Lj_ratio(i,j,ParticleData, dR, dR_check=False):
	"""returns the ratio of hits in Li and Lj when given the data of one particle from the cluster matching function. if dR_check==True, dR is used as a constraint."""
        Li,Lj = 0,0
        for cluster in ParticleData:
		if dR_check and cluster[0][8] >= dR: continue
                if cluster[0][2] == i:
                        Li += 1
                if cluster[0][2] == j:
                        Lj += 1
        if Lj != 0:
                return Li/float(Lj)
        else:
                return 'zero-div'

def Li_Lj_diff(i,j,ParticleData, dR, dR_check=False):
	"""returns the difference of hits in Li and Lj when given the data of one particle from the cluster matching function. if dR_check==True, dR is used as a constraint."""
        Li,Lj = 0,0
        for cluster in ParticleData:
		if dR_check and cluster[0][8] >= dR: continue
                if cluster[0][2] == i:
                        Li += 1
                if cluster[0][2] == j:
                        Lj += 1
        return Li-Lj
	
def Layer_Hist(title,data,dR,minR=0,minPT=0,Save=False):
	"""creates a bar plot of the global amount of matched clusters for each layer and also plots the ratio of consequtive layers. dR here is only used for the title, all cluasters in the file are taken into account. minR (decay vertex R) and minPT can be used for additional constraint."""
        L1,L2,L3,L4, = 0,0,0,0
        if minR>0 or minPT>0:
                add = ' with R>'+str(minR)+' and pT>'+str(minPT)
        else:
                add = ''
        for particle in data:
                for cluster in particle:
                        if cluster[0][7] >= minR and cluster[0][5] >= minPT:
                                if cluster[0][2] == 1: L1 += 1
                                if cluster[0][2] == 2: L2 += 1
                                if cluster[0][2] == 3: L3 += 1
                                if cluster[0][2] == 4: L4 += 1
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
        if Save:
                if minR>0 or minPT>0:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+'R'+str(minR)+'PT'+str(minPT)+title+'.png'
                else:
                        fig2.savefig('HitsPerLayerDR'+str(dR)+title+'.png')
                        print 'saved as HitsPerLayerDR'+str(dR)+title+'.png'
        plt.show()

def Layer_Hist2(title,data,dR,minR=0,minPT=0,Save=False):
	"""creates a bar plot of the global amount of mached clusters for each layer. dR is only used for the tile, all clusters in the file are taken into account. minR (decay vertex R) and minPT can be used for additional constraint."""
        L1,L2,L3,L4, = 0,0,0,0
        if minR>0 or minPT>0:
                add = ' (R>'+str(minR)+', pT>'+str(minPT)+')'
        else:
                add = ''
        for particle in data:
                for cluster in particle:
                        if cluster[0][7] >= minR and cluster[0][5] >= minPT:
                                if cluster[0][2] == 1: L1 += 1
                                if cluster[0][2] == 2: L2 += 1
                                if cluster[0][2] == 3: L3 += 1
                                if cluster[0][2] == 4: L4 += 1
        fig2, ax2 = plt.subplots(1,1,figsize=(5,5))
        fig2.suptitle(title+add+' in dR<'+str(dR))
        ax2.bar([0.5,1.5,2.5,3.5],[L1,L2,L3,L4],align='center')
        ax2.set_ylabel('Clusters')
        ax2.set_xticks([0.5,1.5,2.5,3.5])
        ax2.set_xticklabels(['L1','L2','L3','L4'])
        #plt.tight_layout(pad=2.0,w_pad=0.5,h_pad=0.5)
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

def FillSeparateLayerHist(L1_Hist,L2_Hist,L3_Hist,L4_Hist,title,data,ran,minPT=0):
	"""helper function used in SeparateLayerHist() to fill the histograms"""
        for particle in data:
                L1,L2,L3,L4, = 0,0,0,0
                for n,cluster in enumerate(particle):
			if cluster[0][5] >= minPT:
                        	if cluster[0][2] == 1: L1 += 1
                      		if cluster[0][2] == 2: L2 += 1
                        	if cluster[0][2] == 3: L3 += 1
                        	if cluster[0][2] == 4: L4 += 1
                L1_Hist.Fill(L1)
                L2_Hist.Fill(L2)
                L3_Hist.Fill(L3)
                L4_Hist.Fill(L4)

def SeparateLayerHist(datalist, ran, dR,minPT=0, Save=False):
	"""Creates for each layer and for each dataset given to the datalist as tuple (data, title) a histogram with the amount of clusters in the x-axis"""
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
                FillSeparateLayerHist(Hist_list[0][n],Hist_list[1][n],Hist_list[2][n],Hist_list[3][n],data[1],data[0],ran,minPT=minPT)
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

def ROC_CutBased(title,signal_hist,background_hist,cutregion="above",zerodiv=(0,0),resolution=100):
	"""plots a ROC curve when given a signal and a background histogram. the string cutregion denotes wheter the region 'above' or 'below' should be used for tagging. for a ratio discriminant, the amount of zero-division occurences needs to be given as zerodiv=(nSignal,nBG). If a resolution of 0 is given, it will be set to the amount of integer values inside the histogram range"""
	ran = (signal_hist.GetXaxis().GetXmin(),signal_hist.GetXaxis().GetXmax())
	signal_efficiency = []
	background_efficiency = []
	cuts = np.linspace(ran[0],ran[1],resolution)
	if resolution == 0: cuts=range(int(ran[0]),int(ran[1]))
	bin_ran = (signal_hist.GetXaxis().FindBin(ran[0]),signal_hist.GetXaxis().FindBin(ran[1]))
	if cutregion == "above":
		for cut in cuts[::-1]:
			bin_cut = signal_hist.GetXaxis().FindBin(cut)
			signal_efficiency.append(signal_hist.Integral(bin_cut,bin_ran[1])/(signal_hist.Integral(bin_ran[0],bin_ran[1])+zerodiv[0]))
			background_efficiency.append(1-(background_hist.Integral(bin_cut,bin_ran[1])/(background_hist.Integral(bin_ran[0],bin_ran[1])+zerodiv[1]))) 	
	if cutregion == "below":
		for cut in cuts[::-1]:
			bin_cut = signal_hist.GetXaxis().FindBin(cut)
			signal_efficiency.append(signal_hist.Integral(bin_ran[0],bin_cut)/(signal_hist.Integral(bin_ran[0],bin_ran[1])+zerodiv[0]))
			background_efficiency.append(1-(background_hist.Integral(bin_ran[0],bin_cut)/(background_hist.Integral(bin_ran[0],bin_ran[1])+zerodiv[1]))) 	
	plt.plot(signal_efficiency,background_efficiency,'-',label=title)
	diff = 1
	closest = 0
	for n,BG_eff in enumerate(background_efficiency):
		if abs(BG_eff - 0.9) < diff:
			closest = n
			diff = abs(BG_eff - 0.9)
	print title,"1-BG_efficiency:",background_efficiency[closest], "corresponding to a cut at", cuts[closest]
	diff = 1
	closest = 0
	for n,eff in enumerate(signal_efficiency):
		if abs(eff - 0.6) < diff:
			closest = n
			diff = abs(eff - 0.6)
	print title,"signal_efficiency:",signal_efficiency[closest], "corresponding to a cut at", cuts[closest]
	#fig,ax = plt.subplots()
	#ax.plot(cuts,signal_efficiency,'.', label=r"$\epsilon$ signal")
	#ax.plot(cuts,background_efficiency,'.', label=r"1-$\epsilon$ background")
	#ax.legend(loc=3)
	#ax.set_xlabel("cut")
	#fig.savefig("cutplots/signal_vs_cut_"+title+".png")
	#plt.show()

def Draw_ROC_curves(discriminants,Resolution=0,Title='',ZeroDiv=False,AddCSV=False, pT_hadron=0,pT_jet=0, Save=False):
	"Draws ROC curves for the discriminants given as list of strings ('Li-Lj_dR_X' or 'Li_Lj_dR_X'), which then opens the corresponding root files containing the necessary histograms. an additional title can be given, as well as additional constraints on pT_hadron and pT_jet. If available, the corresponding CSV curve can be added."""
	pT_string = ''
	if pT_hadron>0: pT_string = pT_string+'_pTH_'+str(pT_hadron)
	if pT_jet>0: pT_string = pT_string+'_pTJ_'+str(pT_jet)
	fig1 = plt.figure("2TeV-signal")
	plt.clf()
	fig2 = plt.figure("4TeV-signal")
	plt.clf()
	for discriminant in discriminants:	
		print "working on",discriminant
		if discriminant[2] == '-':
			Resolution = 0
		else:
			Resolution = 50000
		dis_file = rt.TFile("histogram_files/"+discriminant+"histograms.root","READ")
		signal_hist1 = dis_file.Get("2TeV-signal")
		signal_hist2 = dis_file.Get("4TeV-signal")
		background_hist = dis_file.Get("Background")
		if ZeroDiv and discriminant[2]=='_':
			with open("histogram_files/"+discriminant+"_zero_div.pkl",) as f:   
                		zero_div = pickle.load(f)
		else:
			zero_div = {'2TeV-signal':0, '4TeV-signal':0, 'background':0}
		plt.figure("2TeV-signal")
		ROC_CutBased(discriminant,signal_hist1,background_hist,cutregion="above",zerodiv=(zero_div['2TeV-signal'],zero_div['background']),resolution=Resolution)
		plt.figure("4TeV-signal")
		ROC_CutBased(discriminant,signal_hist2,background_hist,cutregion="above",zerodiv=(zero_div['4TeV-signal'],zero_div['background']),resolution=Resolution)
	plt.figure("2TeV-signal")
	if AddCSV:
		print "working on CSV"
		dis_file = rt.TFile("histogram_files/CSV"+pT_string+"histograms.root","READ")
		CSV_hist1 = dis_file.Get("2TeV-signal")
		CSV_hist2 = dis_file.Get("4TeV-signal")
		CSV_background_hist = dis_file.Get("Background")
		ROC_CutBased('CSV',CSV_hist1,CSV_background_hist,cutregion="above",zerodiv=(0,0),resolution=5000)
	plt.title("ROC-curve for 2TeV-signal")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"1-$\epsilon$_background")
	plt.legend(loc=3)
	plt.figure("4TeV-signal")
	if AddCSV:
		ROC_CutBased('CSV',CSV_hist2,CSV_background_hist,cutregion="above",zerodiv=(0,0),resolution=5000)
	plt.title("ROC-curve for 4TeV-signal")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"1-$\epsilon$_background")
	plt.legend(loc=3)
	if Save:
		fig1.savefig("ROC/ROC-curves_2TeV_"+Title+".png")
		print "saved as ROC/ROC-curves_2TeV_"+Title+".png"
		fig2.savefig("ROC/ROC-curves_4TeV_"+Title+".png")
		print "saved as ROC/ROC-curves_4TeV_"+Title+".png"
	plt.show()
	
def tagged_jets_hist(datalist,discriminant, dR, discriminant_cut, CSV_cut, bins, Difference=False, mode="pT_hadron",Save=False):
	"""creates a histogram for each dataset given as list of tuples (data, title, range) of all the jets that were b-tagged by passing a given cut value for CSV and a given discriminant versus a feature given as string to 'mode' (see FEatureDict). The histograms are saved to a root file for further use."""
	title = "tagged_jets_vs_"+mode+"_dR_"+str(dR)
	AllJetsHistlist = []
	CSVHistlist = []
	DiscriminantHistlist = []
	#FeatureDict = {"pT_hadron":5, "pT_jet":10, "decay_vx":7}
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
			AllJetsHistlist[n].Fill(particle[0][0][FeatureDict[mode]])
			if particle[0][0][9] >= CSV_cut: CSVHistlist[n].Fill(particle[0][0][FeatureDict[mode]])
                        if Difference:
                                L = Li_Lj_diff(4,1,particle, dR, dR_check=True)
                        else:
                                L = Li_Lj_ratio(4,1,particle, dR, dR_check=True)
                       	if L == 'zero-div':
                               	continue
			elif L >= discriminant_cut: DiscriminantHistlist[n].Fill(particle[0][0][FeatureDict[mode]])
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
		CSVHistlist[n].Draw("SAME")
		DiscriminantHistlist[n].Draw("SAME")
		legendlist[n].Draw()
		if Save:
			canvaslist[n].SaveAs(title+"_"+data[1]+discriminant+".png")
			Tfilelist.append(rt.TFile("histogram_files/pT_hists/"+title+"_"+data[1]+discriminant+".root","recreate"))
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


def exclusive_tagged_jets_hist(signal_title, data, discriminant, dR, discriminant_cut, CSV_cut,ran, bins, Difference=False, mode="pT_hadron",Save=False):
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
		if particle[0][0][9] >= CSV_cut: 
			CSV_tag = True
                if Difference:
                        L = Li_Lj_diff(4,1,particle, dR, dR_check=True)
                else:
                        L = Li_Lj_ratio(4,1,particle, dR, dR_check=True)
               	if L == 'zero-div':
                       	continue
		elif L >= discriminant_cut: 
			Disc_tag = True
		if Disc_tag and not CSV_tag: Discriminant_and_not_CSV.Fill(particle[0][0][FeatureDict[mode]]) 
		if CSV_tag and not Disc_tag: CSV_and_not_Discriminant.Fill(particle[0][0][FeatureDict[mode]]) 
		if CSV_tag and Disc_tag: Discriminant_and_CSV.Fill(particle[0][0][FeatureDict[mode]]) 
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

def RebinHist(hist,name):
	"""Rebins an efficiency-vs-pT histogram such that the lower statistics at high pT are taken into account"""
	import array
	bins_ = array.array('d',[350.0, 387.5, 425.0, 462.5, 500.0, 537.5, 575.0, 612.5, 650.0, 687.5, 725.0, 762.5, 800.0, 837.5, 875.0, 912.5, 950.0, 987.5, 1025.0, 1100.0, 1200.0, 1400.0, 1600.0, 1800.0, 2100.0, 2500.0, 3000.0])
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

def Make_ROC_histograms(title, data, dR_d, dR_r, pT_cut,Check=False):
	"""obsolete since the Efficient_Cluster_Matcher"""
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
	
	for n,particle in enumerate(data):
		L1_d,L1_r,L4_d,L4_r = 0,0,0,0
		CSV = particle[0][0][9]
		for cluster in particle:
			if cluster[0][8] <= dR_r:
				if cluster[0][2] == 1: L1_r += 1.
				if cluster[0][2] == 4: L4_r += 1.
				if cluster[0][8] <= dR_d:
					if cluster[0][2] == 1: L1_d += 1
					if cluster[0][2] == 4: L4_d += 1
		Diff_hist.Fill(L4_d-L1_d)
		CSV_hist.Fill(CSV)
		if L1_r != 0:
			Ratio_hist.Fill(L4_r/L1_r)
		else:
			ZeroDiv += 1
		if particle[0][0][FeatureDict["pT_hadron"]] >= pT_cut:
			Diff_hist_hadron_pT.Fill(L4_d-L1_d)
			CSV_hist_hadron_pT.Fill(CSV)
			if L1_r != 0:
				Ratio_hist_hadron_pT.Fill(L4_r/L1_r)
			else:
				ZeroDiv_hadron += 1
		if particle[0][0][FeatureDict["pT_jet"]] >= pT_cut:
			Diff_hist_jet_pT.Fill(L4_d-L1_d)
			CSV_hist_jet_pT.Fill(CSV)
			if L1_r != 0:
				Ratio_hist_jet_pT.Fill(L4_r/L1_r)
			else:
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


def Make_ROC_Curves(title,Signal_title,Background_title,diff_ran,ratio_ran,ratio_bins,pT_cut):
	"""Uses the histograms pre-made by efficient_Make_ROC_histograms() to draw ROC-curves in two different representations. The signal and background title should correspond to the ones used in the latter function. The printouts yield which cut corresponds to a 10% mistag rate"""
	diff_bins = diff_ran[1]-diff_ran[0]
	Signal_ZeroDiv = np.loadtxt("histogram_files/{}_ZeroDiv.csv".format(Signal_title),delimiter=',')
	Signal_file = 			rt.TFile("histogram_files/{}_histograms.root".format(Signal_title),"READ")
	Signal_Diff_eff =  		Get_ROC_Efficiencies(Signal_file.Get("L4-L1"),diff_ran,diff_bins,0)
	Signal_Diff_eff_hadron_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4-L1_pTH"),diff_ran,diff_bins,0)
	Signal_Diff_eff_jet_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4-L1_pTJ"),diff_ran,diff_bins,0)
	Signal_Ratio_eff = 		Get_ROC_Efficiencies(Signal_file.Get("L4_L1"),ratio_ran,ratio_bins,Signal_ZeroDiv[0])
	Signal_Ratio_eff_hadron_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4_L1_pTH"),ratio_ran,ratio_bins,Signal_ZeroDiv[1])
	Signal_Ratio_eff_jet_pT = 	Get_ROC_Efficiencies(Signal_file.Get("L4_L1_pTJ"),ratio_ran,ratio_bins,Signal_ZeroDiv[2])
	Signal_CSV_eff = 		Get_ROC_Efficiencies(Signal_file.Get("CSV"),(0,1),ratio_bins,0)
        Signal_CSV_eff_hadron_pT =	Get_ROC_Efficiencies(Signal_file.Get("CSV_pTH"),(0,1),ratio_bins,0)
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
		
	
def ANN_data(signal_data,background_data,dR_d,dR_r):
	"""has become obsolete since the Efficient_Cluster_Matcher() """
	m1 = len(signal_data)
	m2 = len(background_data)
	m = m1+m2
	print "there are {} signal events and {} background events".format(m1,m2)
	csv_file = open("ANN_data.csv","wb")
	writer = csv.writer(csv_file)
	for n,particle in enumerate(signal_data):
		L1_d,L2_d,L3_d,L4_d,L1_r,L2_r,L3_r,L4_r = 0,0,0,0,0,0,0,0
                for cluster in particle:
                        if cluster[0][8] <= dR_d:
                                if cluster[0][2] == 1: L1_d += 1
                                if cluster[0][2] == 2: L2_d += 1
                                if cluster[0][2] == 3: L3_d += 1
                                if cluster[0][2] == 4: L4_d += 1
			elif cluster[0][8] <= dR_r:
                                if cluster[0][2] == 1: L1_r += 1.
                                if cluster[0][2] == 2: L2_r += 1.
                                if cluster[0][2] == 3: L3_r += 1.
                                if cluster[0][2] == 4: L4_r += 1.

		Diff = L4_d - L1_d
		try:
			ratio = L4_r/L1_r
		except ZeroDivisionError:
			ratio = 10*L4_r
		jet_pT = particle[0][0][10]
		CSV = particle[0][0][9]
		writer.writerow([1,CSV,L1_d,L2_d,L3_d,L4_d,Diff,L1_r,L2_r,L3_r,L4_r,ratio,jet_pT])
	
	for n,particle in enumerate(background_data):
		if n >= 2*m1: break
		L1_d,L2_d,L3_d,L4_d,L1_r,L2_r,L3_r,L4_r = 0,0,0,0,0,0,0,0
                for cluster in particle:
                        if cluster[0][8] <= dR_d:
                                if cluster[0][2] == 1: L1_d += 1
                                if cluster[0][2] == 2: L2_d += 1
                                if cluster[0][2] == 3: L3_d += 1
                                if cluster[0][2] == 4: L4_d += 1
			elif cluster[0][8] <= dR_r:
                                if cluster[0][2] == 1: L1_r += 1.
                                if cluster[0][2] == 2: L2_r += 1.
                                if cluster[0][2] == 3: L3_r += 1.
                                if cluster[0][2] == 4: L4_r += 1.

		Diff = L4_d - L1_d
		try:
			ratio = L4_r/L1_r
		except ZeroDivisionError:
			ratio = 10*L4_r
		jet_pT = particle[0][0][10]
		CSV = particle[0][0][9]
		writer.writerow([0,CSV,L1_d,L2_d,L3_d,L4_d,Diff,L1_r,L2_r,L3_r,L4_r,ratio,jet_pT])
	csv_file.close()

def ClusterMatch(title, file_path, dR, MomentumThreshold, JetMode=False, HadronsNotQuarks=False,BG=False, Plot=False, Axes=None, Save=False, dR_dist=False, LightVersion=False, EarlyBreak=0, Continue=False):
        """returns unique ID and coordinates of all pixel clusters that lie inside the dR-cone of a b-particle trajectory; optionally it returns also a 3D-plot
                

        Inputs:
		title:			title used in the filename
                file_path:              path to root file
                dR:                     Delta R region around particle trajectory in which clusters should count as hit
                Momentumthreshold:      momentum threshold for b-particles to be counted
		JetMode:		If set as True, the clusters will be matched with respect to the jet axis instead of the particle trajectory
                HadronsNotQuarks:       if set as True, the algorithm will focus on B-hadrons instead of b-quarks
		BG:			if set as True, only non-signal jets are matched
                Plot:                   if set as True, the function will return a 3D-plot
                Axes:                   pre-initialized 3D-plot axes. Only necessary if Plot==True
                Save:                   if set as True, the function will save the data to a .pkl file for later use
                dR_dist:                if set as True, the function will return a histrogram showing the distribution of delta R between pixel clusters and the corresponding trajectory
		LightVersion:		if set as True, no coordinate data will be stored for clusters in order to save space and computing time
                EarlyBreak:             non zero integer which denotes the number of events after which the algorithm should stop
		Continue:		if set as True, a file with identical parameters will be opened and continued instead of starting from the beginning

        Outputs:
                list of tuples where each contains a tuple with a uniqe identification followed by global cartesian coordinates:((nEvent,nParticle,nLayer,nModule,nCluster,particle_pT,pdgId,decay_vx_R(2D),dR,CSV_b-tag,jet_pT),x,y,z)"""
	
	print "working on file", file_path
        #file = rt.TFile(file_path,"READ")
	file = rt.TFile.Open(file_path)
        #colorstring for plots:
        hsv = plt.get_cmap('hsv')
        color = hsv(np.linspace(0,1.0,12))
        c = 0 #initialize color index
        res = 50 #trajectory resolution

	if JetMode: title=title+'Jet'
	
        # open tree file
        tree = file.Get("demo/tree")
        N = tree.GetEntries()
	print "There are",N,"events in this file."
	if Continue:
		print "opening HitClusterDR"+str(dR)+"on"+title+".pkl"
		with open("HitClusterDR"+str(dR)+"on"+title+".pkl",) as f:
			HitClusters = pickle.load(f)
		N_start = HitClusters[-1][0][0][0]
	else:
        	HitClusters = []
        L1, L2, L3, L4 = 0,0,0,0
        if dR_dist == True:
                histL1 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL2 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL3 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)
                histL4 = rt.TH1D('Clusters(dR) - '+title, 'Clusters(dR) - '+title, 50, 0, dR)

        for i in xrange(N):
		if Continue:
			if i <= N_start: continue
                if i % 100 == 0: print "Working on event " ,i
                if EarlyBreak > 0 and i>=EarlyBreak: break
		if Save and i != 0 and i%10000==0:
			print "saving file - do not abort computation now!" 
                	with open("HitClusterDR"+str(dR)+"on"+title+".pkl", 'w') as f:
                        	pickle.dump(HitClusters, f)
			print "saved as HitClusterDR"+str(dR)+"on"+title+".pkl",
                tree.GetEntry(i)
                for j in range(0,tree.nJets):
                        jVector = rt.TLorentzVector()
                        jVector.SetPtEtaPhiM(tree.jet_pt[j],tree.jet_eta[j],tree.jet_phi[j],tree.jet_mass[j])

			previous_ids = []

                        for k in range(0,tree.nGenParticles):
				if BG:
					pdgCriterion = abs(tree.genParticle_pdgId[k]) != 5# and abs(tree.genParticle_pdgId[k]) < 10
					statusCriterion = tree.genParticle_status[k]== 23
                                elif HadronsNotQuarks == False:
                                        pdgCriterion = abs(tree.genParticle_pdgId[k]) == 5
                                        statusCriterion = tree.genParticle_status[k] == 23
                                else:
                                        pdgCriterion = (abs(tree.genParticle_pdgId[k]) > 500 and abs(tree.genParticle_pdgId[k]) < 600) or (abs(tree.genParticle_pdgId[k]) > 5000 and abs(tree.genParticle_pdgId[k]) < 6000) 
                                        statusCriterion = tree.genParticle_status[k] == 2
                                if statusCriterion and pdgCriterion:
                                        pVector = rt.TLorentzVector()
                                        pVector.SetPtEtaPhiM(tree.genParticle_pt[k],tree.genParticle_eta[k], \
                                                tree.genParticle_phi[k],tree.genParticle_mass[k])
                                        delR = jVector.DeltaR(pVector)
                                        if delR < 0.3 and tree.genParticle_pt[k] > MomentumThreshold: #momentum threshold
                                                if JetMode:
							v_p = normalize(np.array([jVector[0], jVector[1], jVector[2]]))
						else:
							v_p = normalize(np.array([pVector[0], pVector[1], pVector[2]]))
                                                phi = PolarPhi(v_p[0],v_p[1])
                                                theta = Theta(v_p[0],v_p[1],v_p[2])
							
						escape = False          #filter out identical daughters
                                                if len(previous_ids)>0:
                                                        for prid in previous_ids:
                                                                if (abs(abs(prid[0])-abs(tree.genParticle_pdgId[k])) == 2 and abs(prid[0])>100): 
									escape=True
                                                if escape: continue
						
                                                previous_ids.append((tree.genParticle_pdgId[k],delR))

                                                if Plot == True:
                                                        t_max = TrajectoryLength(theta,v_p)
							if JetMode:
								PlotTrajectory((tree.PV_x[0],tree.PV_y[0],tree.PV_z[0]),v_p,Axes,t_max,res,color[c],1,'--')
							else:
                                                        	PlotTrajectory((tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k]),v_p,Axes,t_max,res,color[c],1,'--')

                                                NearClusters = [] #list that will contain for each hit cluster: ((nEvent,nParticle,nModule,nCluster),x,y,z)
                                                for nModule,lenModule in enumerate(tree.nClusters): #finding all clusters inside deltaR<dR
                                                        for nCluster in xrange(0,lenModule):
								if JetMode:
									ClusterTheta,ClusterPhi = ShiftedThetaPhi(tree.cluster_globalx[nModule][nCluster],\
                                                                        	tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster],\
                                                                        	tree.PV_x[0],tree.PV_y[0],tree.PV_z[0])

								else:
									ClusterTheta,ClusterPhi = ShiftedThetaPhi(tree.cluster_globalx[nModule][nCluster],\
                                                                        	tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster],\
                                                                        	tree.genParticle_vx_x[k],tree.genParticle_vx_y[k],tree.genParticle_vx_z[k])
                                                                DR = DeltaR(theta,ClusterTheta,phi,ClusterPhi)
                                                                if DR<dR:
									if LightVersion:
										NearClusters.append(((i,k,tree.detUnit_layer[nModule],nModule,nCluster,tree.genParticle_pt[k],tree.genParticle_pdgId[k],np.sqrt(tree.genParticle_decayvx_x[k]**2+tree.genParticle_decayvx_y[k]**2),DR,tree.jet_bTag[j],tree.jet_pt[j]),''))
 									else:	
                                                                        	NearClusters.append(((i,k,tree.detUnit_layer[nModule],nModule,nCluster,tree.genParticle_pt[k],tree.genParticle_pdgId[k],np.sqrt(tree.genParticle_decayvx_x[k]**2+tree.genParticle_decayvx_y[k]**2),DR,tree.jet_bTag[j],tree.jet_pt[j]),tree.cluster_globalx[nModule][nCluster],tree.cluster_globaly[nModule][nCluster],tree.cluster_globalz[nModule][nCluster]))

                                                                        if dR_dist == True:
                                                                                if tree.detUnit_layer[nModule] == 1:
                                                                                        histL1.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 2:
                                                                                        histL2.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 3:
                                                                                        histL3.Fill(DR)
                                                                                elif tree.detUnit_layer[nModule] == 4:
                                                                                        histL4.Fill(DR)
                                                                        if LayerHist == True:
                                                                                if tree.detUnit_layer[nModule] == 1:
                                                                                        L1 += 1
                                                                                elif tree.detUnit_layer[nModule] == 2:
                                                                                        L2 += 2
                                                                                elif tree.detUnit_layer[nModule] == 3:
                                                                                        L3 += 3
                                                                                elif tree.detUnit_layer[nModule] == 4:
                                                                                        L4 += 4
                                                
                                                if  NearClusters == []:
                                                        break 
                                                else:
                                                        HitClusters.append(NearClusters) #summarizing the hit cluster ID for every event
                                                        if Plot == True:
								X,Y,Z=[],[],[]
								for entry in NearClusters:
                                                                	X.append(entry[1])
                                                                	Y.append(entry[2])
                                                                	Z.append(entry[3])

                                                                Axes.scatter(X,Y,Z,c=color[c],s=9,linewidths=0.1) #plots all the hit clusters
                                                                if c != len(color)-1:
                                                                        c += 1
                                                                else:
                                                                        c = 0
        if BG: 
		particleString = 'background-particles'
	elif HadronsNotQuarks == False:
                particleString = 'b-quarks'
        else:
                particleString = 'B-hadrons'
        print "Total Number of high pt "+particleString+": ", len(HitClusters)
        nHitClusters = 0
        for entry in HitClusters:
                nHitClusters += len(entry)
        print "Total number of clusters hit: ", nHitClusters
        if Plot == True:
                plt.savefig("Trajectories_Clusters.png")
                plt.show()
        if Save == True:
		print "saving file - do not abort computation now!" 
                with open("HitClusterDR"+str(dR)+"on"+title+".pkl", 'w') as f:
                        pickle.dump(HitClusters, f)
		print "Saved as HitClusterDR"+str(dR)+"on"+title+".pkl"

        if dR_dist == True:
                c1 = rt.TCanvas('c1','c1',600,600)
                c1.SetTitle('dR distribution')
		rt.gStyle.SetOptStat(0)
                histL1.GetXaxis().SetTitle('dR')
                histL1.GetYaxis().SetTitle('[a.u.]')
                histL1.GetYaxis().SetTitleOffset(1.5)
                histL1.SetLineColor(1)
                histL2.SetLineColor(2)
                histL3.SetLineColor(3)
                histL4.SetLineColor(4)
                histL1.DrawNormalized()
                histL2.DrawNormalized('SAME')
                histL3.DrawNormalized('SAME')
                histL4.DrawNormalized('SAME')
                if BG:
			#l1 = rt.TLegend(0.9,0.1,0.65,0.25)
			l1 = rt.TLegend(0.9,0.9,0.65,0.75)
		else:
			l1 = rt.TLegend(0.9,0.9,0.65,0.75)
                l1.AddEntry(histL1,'Layer 1')
                l1.AddEntry(histL2,'Layer 2')
                l1.AddEntry(histL3,'Layer 3')
                l1.AddEntry(histL4,'Layer 4')
                l1.Draw()
                c1.SaveAs('dR_dist/DeltaR-dist'+title+'2.png')
		#print 'saved as dR_dist/DeltaR-dist'+title+'.png'

	return HitClusters

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
                creates numpy array where each row correspond to a jet with the attributes in the following order: nEvent, CSV-tag, hadron_pT, jet_pT, decay_vx, Li_0.04, Li_0.06, Li_0.08, Li_0.1, Li_0.16 """

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
		HitClusters = np.load("MatchedClusters_{}.npy".format(title))
		N_start = int(HitClusters[-1,0]+1)
		print "successfully loaded MatchedClusters_{}.npy and continue on event {}".format(title,N_start) 
	else:
        	HitClusters = np.ndarray(shape=(0,25))
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
			np.save("MatchedClusters_{}.npy".format(title),HitClusters)	
			print "saved as MatchedClusters_{}.npy".format(title)
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
									
						HitClusters = np.vstack((HitClusters,np.concatenate((JetFeatures,HitsPerLayer))))

        print "Total Number of matched high pt jets:",HitClusters.shape[0]
	print "saving file - do not abort computation now!" 
	np.save("MatchedClusters_{}.npy".format(title),HitClusters)
	print "Saved as Matched_Clusters_{}.npy".format(title)
	t2 = time.time()
	print "total processing time: {}s".format(t2-t1)
	if Protocol:
		writer.writerow([title,N,N])
		csv_file.close()


if __name__ == '__main__':

	#select global parameters

        #dR = 0.16 #DeltaR threshold for counting clusters
        MomentumThreshold = 350


	'''
	#initialize 3D-plot

        with open("Grid.pkl",) as f:    #open coordinates of DetUnits for visual reference
                Grid = pickle.load(f)
	ax = Initialize3DPlot('Particle_Trajectories', 'x', 'y', 'z', grid=Grid)
	'''
	'''
	print "HitClusterDR0.1onDummy_4TeV-SignalJet.pkl" 	
	with open("HitClusterDR0.1onDummy_4TeV-SignalJet.pkl",) as f:   
                Signal = pickle.load(f)

	print "HitClusterDR0.16onBackgroundJet.pkl"
	with open("HitClusterDR0.16onBackgroundJet.pkl",) as f: 
               Background = pickle.load(f)
	'''

	
	'''
	print "opening file HitClusterDR0.16on2TeV-SignalJet.pkl"
	with open("HitClusterDR0.16on2TeV-SignalJet.pkl",) as f:   
                Signal1 = pickle.load(f)
	
	print "opening file HitClusterDR0.16on4TeV-SignalJet.pkl" 	
	with open("HitClusterDR0.16on4TeV-SignalJet.pkl",) as f:   
                Signal2 = pickle.load(f)
	'''
	#print "opening file HitClusterDR0.1onBG1Jet.pkl"
	#with open("HitClusterDR0.1onBG1Jet.pkl",) as f: 
        #       Background = pickle.load(f)
	
	
	#ANN_data(Signal2,Background,0.04,0.1)
	#Make_ROC_histograms("4TeV-Signal", Signal2, 0.04, 0.01, 1200,Check=False)
	'''
	Background = []
	for n in range(1,4):
		print "opening file HitClusterDR0.1onBG{}Jet.pkl".format(n)
		with open("HitClusterDR0.1onBG{}Jet.pkl".format(n),) as f:
			Background += pickle.load(f)
	'''
		
	'''
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

	
	
	#Correlation_Hist2D("jet_pT_vs_decayvx_R", Signal2, "pT_jet","decay_vx", (0,2500), (0,60), 50, 50,Profiling = False,Save=True)
	#Correlation_Hist2D("jet_pT_vs_decayvx_R_Profiling", Signal2, "pT_jet","decay_vx", (0,2500), (0,60), 50, 50,Profiling = True,Save=True)
	
	#Correlation_Hist2D("jet_pT_vs_hadron_pT", Signal2, "pT_jet","pT_hadron", (0,2500), (0,2500), 50, 50,Profiling = False,Save=True)
	#Correlation_Hist2D("jet_pT_vs_hadron_pT_Profiling", Signal2, "pT_jet","pT_hadron", (0,2500), (0,2500), 50, 50,Profiling = True,Save=True)

	'''
	'''
	#Global HitCluster count per layer

	Layer_Hist2('2TeV',Signal1,dR=dR,minR=0, minPT=0, Save=True)
	Layer_Hist2('4TeV',Signal2,dR=dR, minR=0, minPT=0, Save=True)
	Layer_Hist2('Background',Background,dR=dR, minR=0, minPT=0, Save=True)

	Layer_Hist2('2TeV',Signal1,dR=dR,minR=4, minPT=1000, Save=True)
	Layer_Hist2('4TeV',Signal2,dR=dR, minR=4, minPT=1000, Save=True)
	Layer_Hist2('Background',Background,dR=dR, minR=4, minPT=1000, Save=True)
	'''
        
	#Separate HitCluster per Layer Histograms

	#SeparateLayerHist([(Signal1,'2Tev-signal'),(Signal2,'4Tev-signal'),(Background,'Background')], (0,30), dR, minPT=0, Save=True)

        
	#select file paths	

	SignalFile1 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M2000_GENSIMDIGIRECO_v2.root"
	#SignalFile2 = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/ZprimeBBbar_M4000_GENSIMDIGIRECO_v2.root"
	#Background = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8.root"
	#Background = "/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/QCD_noPU_v2.root"
	
	#MakeGrid(SignalFile1,EarlyBreak=2000)
	
	#additional Background files

	Additional_Background_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_wPVs/180502_130824/0000/flatTuple_{}.root'
	
	#pre-process data

	#Bdata1 =  ClusterMatch('2TeV-Signal', SignalFile1, dR, MomentumThresholdB,JetMode=True, HadronsNotQuarks=True, Plot=False, Axes=None, Save=False, dR_dist = False, LayerHist=False, LightVersion=True, EarlyBreak=1)
	#Bdata2 =  ClusterMatch('4TeV-Signal', SignalFile2, dR, MomentumThresholdB,JetMode=True, HadronsNotQuarks=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False,LightVersion=True, EarlyBreak=0)
	#Backgrounddata = ClusterMatch('Background', Background, dR, MomentumThresholdBackground,JetMode=True, HadronsNotQuarks=True, BG=True, Plot=False, Axes=None, Save=True, dR_dist = False, LayerHist=False,LightVersion=True, EarlyBreak=180000)
	#Efficient_Cluster_Matcher("4TeV-Signal", SignalFile2, MomentumThreshold, BG=False, EarlyBreak=0, Continue=True)
	#Efficient_Cluster_Matcher("Background", Background, MomentumThreshold, BG=True, EarlyBreak=0, Continue=True,Protocol=True)
	Efficient_Cluster_Matcher("2TeV-Signal", SignalFile1, 350, BG=False, EarlyBreak=0, Continue=False, Protocol=False)




	'''	
	#need to get GRID-permission first!!!
	for n in range(23,30):
		try:
			Efficient_Cluster_Matcher("BG"+str(n), Additional_Background_String.format(n), MomentumThreshold, BG=True, EarlyBreak=0, Continue=True,Protocol=True)
		except:
			continue
	for n in range(30,35):
		try:
			Efficient_Cluster_Matcher("BG"+str(n), Additional_Background_String.format(n), MomentumThreshold, BG=True, EarlyBreak=0, Continue=False,Protocol=True)
		except:
			continue

	csv_file = open("Protocol.csv","wb")
       	writer = csv.writer(csv_file)
	writer.writerow([0,0,0])
	csv_file.close()
	'''
	'''
	#HugeBackground = ClusterMatch('Background', Additional_Background[0], dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=True, Axes=ax, Save=False, dR_dist = False, LayerHist=False, EarlyBreak=500)
	

	'''
	
	#Plot CSV-histograms

	#CSV_Hist([(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, pT_hadron=1000, Save=True)
	#CSV_Hist([(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, pT_jet=1000, Save=True)
	#Draw_ROC_curves(['CSV'], Resolution = 50000,Title='CSV',Save=True)
	
	'''
	
	#Plot 2D-Histograms for every combination of Li,Lj and for every available sample
	
	Li_Lj_Hist2D('2TeV-Signal',1,2,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',1,3,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',1,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',2,3,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',2,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('2TeV-Signal',3,4,Signal1,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,2,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,3,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',1,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',2,3,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',2,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('4TeV-Signal',3,4,Signal2,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,2,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,3,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',1,4,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',2,3,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',2,4,Background,(0,35),dR,Save=True)
	Li_Lj_Hist2D('background',3,4,Background,(0,35),dR,Save=True)
	'''	
	'''
	#Plot 1D-Histogram for different combinations of Li,Lj
	
	Li_Lj_Hist1D(2, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 3, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L2-L1_dR_0.02','L3-L1_dR_0.02','L4-L1_dR_0.02','L3-L2_dR_0.02','L4-L2_dR_0.02','L4-L3_dR_0.02']
	Draw_ROC_curves(hist_files, Resolution = 0,Title='Li-Lj_dR_0.02_Jet',AddCSV=True,Save=False)
	
	Li_Lj_Hist1D(2, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 3, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L2-L1_dR_0.04','L3-L1_dR_0.04','L4-L1_dR_0.04','L3-L2_dR_0.04','L4-L2_dR_0.04','L4-L3_dR_0.04']
	Draw_ROC_curves(hist_files, Resolution = 0,Title='Li-Lj_dR_0.04_JetCSV',AddCSV=True,Save=False)
	
	Li_Lj_Hist1D(2, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 3, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L2_L1_dR_0.04','L3_L1_dR_0.04','L4_L1_dR_0.04','L3_L2_dR_0.04','L4_L2_dR_0.04','L4_L3_dR_0.04']
	Draw_ROC_curves(hist_files, Resolution = 50000,Title='Li_Lj_dR_0.04_JetCSV',ZeroDiv=True,AddCSV=True,Save=False)
	
	Li_Lj_Hist1D(2, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(3, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 2, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 3, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L2_L1_dR_0.08','L3_L1_dR_0.08','L4_L1_dR_0.08','L3_L2_dR_0.08','L4_L2_dR_0.08','L4_L3_dR_0.08']
	Draw_ROC_curves(hist_files, Resolution = 50000,Title='Li_Lj_dR_0.08_Jet2',ZeroDiv=True,AddCSV=True,Save=False)
	
	
	
	#Plot 1D-Histogram for L4_L1 and different dR
	
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.01,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.16,Difference=True, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L4-L1_dR_0.01', 'L4-L1_dR_0.02', 'L4-L1_dR_0.04', 'L4-L1_dR_0.08', 'L4-L1_dR_0.16'] 
	Draw_ROC_curves(hist_files, Resolution = 0,Title='L4-L1_Jet',AddCSV=True,Save=False)
	
	
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.01,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.02,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.08,Difference=False, Abs=False, dR_check=True, Save=True)
	
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.1,Difference=False, Abs=False, dR_check=True, Save=True)
	
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.16,Difference=False, Abs=False, dR_check=True, Save=True)
	
	hist_files = ['L4_L1_dR_0.01', 'L4_L1_dR_0.02', 'L4_L1_dR_0.04', 'L4_L1_dR_0.08','L4_L1_dR_0.1','L4_L1_dR_0.16'] 
	Draw_ROC_curves(hist_files, Resolution = 50000,Title='L4_L1_Jet2',ZeroDiv=True,AddCSV=True,Save=False)
	
	
	#Best discriminants

	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,pT_hadron=1000, Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.1,pT_hadron=1000, Difference=False, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.04,pT_jet=1000, Difference=True, Abs=False, dR_check=True, Save=True)
	Li_Lj_Hist1D(4, 1, [(Signal1,'2TeV-signal'),(Signal2,'4TeV-signal')], Background, 40, (0,6),dR=0.1,pT_jet=1000, Difference=False, Abs=False, dR_check=True, Save=True)

	'''
	'''	
	hist_files = ['L4-L1_dR_0.04_pTH_1200','L4-L1_dR_0.04_pTJ_1200','L4_L1_dR_0.1_pTH_1200','L4_L1_dR_0.1_pTJ_1200','CSV_pTH_1200','CSV_pTJ_1200']
	hist_files += ['L4-L1_dR_0.04','L4_L1_dR_0.1','CSV']
	#hist_files = ['L4-L1_dR_0.04','L4_L1_dR_0.1','CSV', 'L4-L1_dR_0.04_pTJ_1000','L4_L1_dR_0.1_pTJ_1000','CSV_pTJ_1000']
	Draw_ROC_curves(hist_files, Title='Best_Discriminants',ZeroDiv=True,AddCSV=False,Save=True)
	'''

	
	'''
	#manually print out all the clusters hit by background particles

	BG1 = ClusterMatch(,'BG1',Background_file, dR, MomentumThresholdBackground, HadronsNotQuarks=True, BG=True, Plot=True, Axes=ax, Save=False, dR_dist = False, LayerHist=False, EarlyBreak=500)
	for particle in BG1:
		print "particle id =", particle[0][0][6]
		print "particle pT", particle[0][0][5]
		L1, L2, L3, L4, LB = 0,0,0,0,0
		for cluster in particle:
			if cluster[0][2] == 1: L1+=1
			if cluster[0][2] == 2: L2+=1
			if cluster[0][2] == 3: L3+=1
			if cluster[0][2] == 4: L4+=1
			if cluster[0][2] == 0: LB +=1
		print "Hits in L1 =",L1,", L2 =",L2,", L3 =",L3,"L4 =",L4,"Barrel =",LB
	'''
	
	#L2_L1 = lambda ParticleData : Li_Lj(2,1,ParticleData)
	#L3_L1 = lambda ParticleData : Li_Lj(3,1,ParticleData)
        #L4_L1 = lambda ParticleData : Li_Lj(4,1,ParticleData)
	#L3_L2 = lambda ParticleData : Li_Lj(3,2,ParticleData)
	#L4_L2 = lambda ParticleData : Li_Lj(4,2,ParticleData)
	#L4_L3 = lambda ParticleData : Li_Lj(4,3,ParticleData)
	#DiscriminatorHist('L2_L1',L2_L1,Bdata,Backgrounddata,40,(0,6),'L2/L1')
	#DiscriminatorHist('L3_L1',L3_L1,Bdata,Backgrounddata,40,(0,6),'L3/L1')
	#DiscriminatorHist('L4_L1',L4_L1,Bdata,Backgrounddata,40,(0,6),'L4/L1')
	#DiscriminatorHist('L3_L2',L3_L2,Bdata,Backgrounddata,40,(0,6),'L3/L2')
	#DiscriminatorHist('L4_L2',L4_L2,Bdata,Backgrounddata,40,(0,6),'L4/L2')
	#DiscriminatorHist('L4_L3',L4_L3,Bdata,Backgrounddata,40,(0,6),'L4/L3')
	
	
	'''
	#search for ideal status requirement on background particles
	
	tree = Background_file.Get("demo/tree")
        nBE = tree.GetEntries()
	#print "number of B events: ", nBE
	#print "Total Number of high pt B-hadrons ", len(Bdata)
        #nB = 0
        for event in range(nBE):
		if event > 500: break
        	tree.GetEntry(event) 
		for particle in range(0,tree.nGenParticles):
			if not (abs(tree.genParticle_pdgId[particle]) < 7 or abs(tree.genParticle_pdgId[particle])==21 or tree.genParticle_status[particle]>10):
				print "particle Id =", tree.genParticle_pdgId[particle]
				print "particle status =", tree.genParticle_status[particle]
	#print "Total number of clusters hit by B: ", nB
	'''
	


	#Histogram of all statuses

	#statusHist = rt.TH1D('status','status',110,0,109)
	#tree = Background_file.Get("demo/tree")
        #nBGE = tree.GetEntries()
	#for i in xrange(nBGE):
	#	tree.GetEntry(i)
	#	if i%100 == 0: print "working on event", i
	#	for particle in range(0,tree.nGenParticles):
        #               statusHist.Fill(tree.genParticle_status[particle])
	#canvas = rt.TCanvas('canvas','canvas',600,600)
        #canvas.SetTitle("status")
        #statusHist.GetXaxis().SetTitle("status")
        #statusHist.GetYaxis().SetTitle('# particles')
        #statusHist.GetYaxis().SetTitleOffset(1.5)
        #statusHist.SetLineColor(2)
        #statusHist.Draw()
        #canvas.SaveAs('BGstatus.png')
	#print "number of background events: ", nBGE
	#print "Total Number of Background particles ", len(Backgrounddata)
        #nBG = 0
        #for entry in Backgrounddata:
        #        nBG += len(entry)
        #print "Total number of clusters hit by background: ", nBG






