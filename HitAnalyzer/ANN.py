import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
import keras as kr
import pickle
import csv
from sklearn.metrics import roc_curve

print "packages imported"


#different models

def build_model():
	model = kr.models.Sequential()

	model.add(kr.layers.Dense(150, activation = 'relu', input_shape=(20,)))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(1, activation ='sigmoid'))
	
	model.compile(optimizer='adagrad',loss='binary_crossentropy',metrics=['accuracy'])
	return model

def build_conv_model():
	model = kr.models.Sequential()

	model.add(kr.layers.Conv2D(64,kernel_size=(5,2),strides=1, border_mode='valid',activation = 'relu', input_shape=(5,4,1)))
	model.add(kr.layers.Flatten())
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(150, activation = 'relu'))
	model.add(kr.layers.Dropout(0.5))
	
	model.add(kr.layers.Dense(1, activation ='sigmoid'))
	
	model.compile(optimizer='adagrad',loss='binary_crossentropy',metrics=['accuracy'])
	return model

def build_functional_model():
	Li_inputs = kr.layers.Input(shape=(5,4,1))
	Li_Lj_inputs = kr.layers.Input(shape=(5,))
	
	Li_branch =  kr.layers.Conv2D(32,kernel_size=(5,2),strides=1, border_mode='valid',activation = 'relu', input_shape=(5,4,1))(Li_inputs)
	Li_branch = kr.layers.Flatten()(Li_branch)

	Li_Lj_branch = kr.layers.Dense(150, activation='relu')(Li_Lj_inputs)
	Li_Lj_branch = kr.layers.Dropout(0.5)(Li_Lj_branch)

	x = kr.layers.concatenate([Li_branch, Li_Lj_branch],axis=1)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)

	output = kr.layers.Dense(1, activation ='sigmoid')(x)	 

	model = kr.models.Model(inputs=[Li_inputs, Li_Lj_inputs], outputs=output)
	model.compile(optimizer='adagrad',loss='binary_crossentropy',metrics=['accuracy'])
	return model

def build_complex_functional_model(Additional_Input=False):
	Li_inputs = kr.layers.Input(shape=(5,4,1))
	Li_Lj_inputs = kr.layers.Input(shape=(5,))
	
	conv_branch =  kr.layers.Conv2D(64,kernel_size=(5,2),strides=1, border_mode='valid',activation = 'relu', input_shape=(5,4,1))(Li_inputs)
	conv_branch = kr.layers.Flatten()(conv_branch)
	
	Li_branch = kr.layers.Flatten(input_shape=(5,4,1))(Li_inputs)
	Li_branch = kr.layers.Dense(150, activation='relu')(Li_branch)

	Li_Lj_branch = kr.layers.Dense(150, activation='relu')(Li_Lj_inputs)
	Li_Lj_branch = kr.layers.Dropout(0.5)(Li_Lj_branch)

	if Additional_Input:
		add_input = kr.layers.Input(shape=(1,))
		x = kr.layers.concatenate([conv_branch, Li_branch, Li_Lj_branch, add_input],axis=1)
	else:
		x = kr.layers.concatenate([conv_branch, Li_branch, Li_Lj_branch],axis=1)
	
	x = kr.layers.Dense(200, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)
	x = kr.layers.Dense(150, activation='relu')(x)
	x = kr.layers.Dropout(0.5)(x)

	output = kr.layers.Dense(1, activation ='sigmoid')(x)	 
	if Additional_Input:
		model = kr.models.Model(inputs=[Li_inputs, Li_Lj_inputs, add_input], outputs=output)
	else:
		model = kr.models.Model(inputs=[Li_inputs, Li_Lj_inputs], outputs=output)
	model.compile(optimizer='adagrad',loss='binary_crossentropy',metrics=['accuracy'])
	return model


#data processing functions

def New_Data(title): 
	from FinalClusterMatcher import load_data
	Signal_data = load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
	Background_data = load_data('BG_PU','matched_clusters/BG_PU/',499)
	n = Signal_data.shape[0]
	if Background_data.shape[0] >1.5*n:
		Background_data = Background_data[:int(1.5*n),:]
	Signal_data[:,0]=1
	Background_data[:,0]=0
	Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
	np.random.shuffle(Data_sample)
	m = Data_sample.shape[0]
	r = int(0.2*m)
	np.save("ANN_data/test_x_{}.npy".format(title),Data_sample[:r,3:])
	np.save("ANN_data/test_y_{}.npy".format(title),Data_sample[:r,0].astype(int))
	np.save("ANN_data/test_CSV_{}.npy".format(title),Data_sample[:r,1])
	np.save("ANN_data/test_pT_{}.npy".format(title),Data_sample[:r,2])
	np.save("ANN_data/train_x_{}.npy".format(title),Data_sample[r:,3:])
	np.save("ANN_data/train_y_{}.npy".format(title),Data_sample[r:,0].astype(int))
	np.save("ANN_data/train_pT_{}.npy".format(title),Data_sample[r:,2])


def New_Data_4_Submit(): 
	from FinalClusterMatcher import load_data
	job_dict = {0:'noPU_both', 1:'noPU_both_withPT', 2:'noPU_4TeV', 3:'noPU_4TeV_withPT', 4:'withPU_both', 5:'withPU_both_withPT', 6:'withPU_both_withPV', 7:'withPU_4TeV', 8:'withPU_4TeV_withPT', 9:'withPU_4TeV_withPV'}

	Signal_4TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_4TeV-Signal.npy')
        Signal_2TeV_noPU = np.load('matched_clusters/Signal_noPU/MatchedClusters_2TeV-Signal.npy')
        Signal_both_noPU = np.vstack((Signal_4TeV_noPU,Signal_2TeV_noPU))
        Background_noPU = load_data('BG','matched_clusters/BG_noPU/',31, old_version=True)
        Signal_4TeV_PU = load_data('4TeV-Signal_PU','matched_clusters/Signal_PU/',15)
        Signal_2TeV_PU = load_data('2TeV-Signal_PU','matched_clusters/Signal_PU/',19)
        Signal_both_PU = np.vstack((Signal_4TeV_PU,Signal_2TeV_PU))
        Background_PU = load_data('BG_PU','matched_clusters/BG_PU/',499)

	#noPU_4TeV
	Signal_data = Signal_4TeV_noPU
	Background_data = Background_noPU

	n = Signal_data.shape[0]
	if Background_data.shape[0] >1.3*n:
		Background_data = Background_data[:int(1.3*n),:]
	Signal_data[:,0]=1
	Background_data[:,0]=0
	Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
	np.random.shuffle(Data_sample)
	m = Data_sample.shape[0]
	r = int(0.2*m)
	for title in [job_dict[2],job_dict[3]]:
		np.save("Submitted_Models/data/{}/test_x.npy".format(title),Data_sample[:r,3:])
		np.save("Submitted_Models/data/{}/test_y.npy".format(title),Data_sample[:r,0].astype(int))
		np.save("Submitted_Models/data/{}/test_CSV.npy".format(title),Data_sample[:r,1])
		np.save("Submitted_Models/data/{}/train_x.npy".format(title),Data_sample[r:,3:])
		np.save("Submitted_Models/data/{}/train_y.npy".format(title),Data_sample[r:,0].astype(int))
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[3]),Data_sample[:r,2])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[3]),Data_sample[r:,2])

	#noPU_both

	Signal_data = Signal_both_noPU
	Background_data = Background_noPU

	n = Signal_data.shape[0]
	if Background_data.shape[0] >1.3*n:
		Background_data = Background_data[:int(1.3*n),:]
	Signal_data[:,0]=1
	Background_data[:,0]=0
	Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
	np.random.shuffle(Data_sample)
	m = Data_sample.shape[0]
	r = int(0.2*m)
	for title in [job_dict[0],job_dict[1]]:
		np.save("Submitted_Models/data/{}/test_x.npy".format(title),Data_sample[:r,3:])
		np.save("Submitted_Models/data/{}/test_y.npy".format(title),Data_sample[:r,0].astype(int))
		np.save("Submitted_Models/data/{}/test_CSV.npy".format(title),Data_sample[:r,1])
		np.save("Submitted_Models/data/{}/train_x.npy".format(title),Data_sample[r:,3:])
		np.save("Submitted_Models/data/{}/train_y.npy".format(title),Data_sample[r:,0].astype(int))
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[1]),Data_sample[:r,2])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[1]),Data_sample[r:,2])

	#withPU_4TeV

	Signal_data = Signal_4TeV_PU
	Background_data = Background_PU

	n = Signal_data.shape[0]
	if Background_data.shape[0] >1.3*n:
		Background_data = Background_data[:int(1.3*n),:]
	Signal_data[:,0]=1
	Background_data[:,0]=0
	Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
	np.random.shuffle(Data_sample)
	m = Data_sample.shape[0]
	r = int(0.2*m)
	for title in [job_dict[7],job_dict[8], job_dict[9]]:
		np.save("Submitted_Models/data/{}/test_x.npy".format(title),Data_sample[:r,3:])
		np.save("Submitted_Models/data/{}/test_y.npy".format(title),Data_sample[:r,0].astype(int))
		np.save("Submitted_Models/data/{}/test_CSV.npy".format(title),Data_sample[:r,1])
		np.save("Submitted_Models/data/{}/train_x.npy".format(title),Data_sample[r:,3:])
		np.save("Submitted_Models/data/{}/train_y.npy".format(title),Data_sample[r:,0].astype(int))
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[8]),Data_sample[:r,2])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[8]),Data_sample[r:,2])
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[9]),Data_sample[:r,-1])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[9]),Data_sample[r:,-1])

	#withPU_both

	Signal_data = Signal_both_PU
	Background_data = Background_PU

	n = Signal_data.shape[0]
	if Background_data.shape[0] >1.3*n:
		Background_data = Background_data[:int(1.3*n),:]
	Signal_data[:,0]=1
	Background_data[:,0]=0
	Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
	np.random.shuffle(Data_sample)
	m = Data_sample.shape[0]
	r = int(0.2*m)
	for title in [job_dict[4],job_dict[5], job_dict[6]]:
		np.save("Submitted_Models/data/{}/test_x.npy".format(title),Data_sample[:r,3:])
		np.save("Submitted_Models/data/{}/test_y.npy".format(title),Data_sample[:r,0].astype(int))
		np.save("Submitted_Models/data/{}/test_CSV.npy".format(title),Data_sample[:r,1])
		np.save("Submitted_Models/data/{}/train_x.npy".format(title),Data_sample[r:,3:])
		np.save("Submitted_Models/data/{}/train_y.npy".format(title),Data_sample[r:,0].astype(int))
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[5]),Data_sample[:r,2])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[5]),Data_sample[r:,2])
	np.save("Submitted_Models/data/{}/test_feature.npy".format(job_dict[6]),Data_sample[:r,-1])
	np.save("Submitted_Models/data/{}/train_feature.npy".format(job_dict[6]),Data_sample[r:,-1])

	'''
	for i in range(10):
		title = job_dict[i]
		if i<4:
			Background_data = Background_noPU
			if i==0 or i==1:
				Signal_data = Signal_both_noPU
			else:
				Signal_data = Signal_4TeV_noPU
		else:
			Background_data = Background_PU
			if i==4 or i==5 or i==6:
				Signal_data = Signal_both_PU
			else:
				Signal_data = Signal_4TeV_PU

		n = Signal_data.shape[0]
		if Background_data.shape[0] >1.3*n:
			Background_data = Background_data[:int(1.3*n),:]
		Signal_data[:,0]=1
		Background_data[:,0]=0
		Data_sample = np.delete(np.vstack((Signal_data,Background_data)),[2,4],axis=1)
		np.random.shuffle(Data_sample)
		m = Data_sample.shape[0]
		r = int(0.2*m)
		np.save("Submitted_Models/data/{}/test_x.npy".format(title),Data_sample[:r,3:])
		np.save("Submitted_Models/data/{}/test_y.npy".format(title),Data_sample[:r,0].astype(int))
		np.save("Submitted_Models/data/{}/test_CSV.npy".format(title),Data_sample[:r,1])
		np.save("Submitted_Models/data/{}/train_x.npy".format(title),Data_sample[r:,3:])
		np.save("Submitted_Models/data/{}/train_y.npy".format(title),Data_sample[r:,0].astype(int))
		
		if i==1 or i==3 or i==5 or i==8:
			np.save("Submitted_Models/data/{}/test_feature.npy".format(title),Data_sample[:r,2])
			np.save("Submitted_Models/data/{}/train_feature.npy".format(title),Data_sample[r:,2])
		elif i==6 or i==9:
			np.save("Submitted_Models/data/{}/test_feature.npy".format(title),Data_sample[:r,-1])
			np.save("Submitted_Models/data/{}/train_feature.npy".format(title),Data_sample[r:,-1])
	'''		

def discriminants(x_data):
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
        x_data_Li_Lj= discriminants(x_data)
        x_data_Li = np.reshape(x_data[:,:20].flatten(),(-1,5,4,1))
        return [x_data_Li, x_data_Li_Lj]

#evaluation functions

def Compare_ANNs():
	fpr_csv, tpr_csv, thresholds_csv = roc_curve(test_y,test_CSV)
	tpr1,fpr1 = np.loadtxt("ANN_data/Only_Li_d.csv",delimiter=',')
	tpr2,fpr2 = np.loadtxt("ANN_data/Only_Li_r.csv",delimiter=',')
	tpr3,fpr3 = np.loadtxt("ANN_data/VariableCone.csv",delimiter=',')
	tpr4,fpr4 = np.loadtxt("ANN_data/Li_d_Diff.csv",delimiter=',')
	tpr5,fpr5 = np.loadtxt("ANN_data/Li_r_Ratio.csv",delimiter=',')
	tpr7,fpr7 = np.loadtxt("ANN_data/Li_d_Diff_Li_r_Ratio.csv",delimiter=',')
	tpr8,fpr8 = np.loadtxt("ANN_data/5Cone.csv",delimiter=',')
	tpr9,fpr9 = np.loadtxt("ANN_data/5Cone_conv.csv",delimiter=',')
	
	plt.figure("ROC")
	plt.clf()
	plt.semilogy(tpr1,fpr1,'r-',label='Li dR<0.04')
	plt.semilogy(tpr2,fpr2,'r:',label='Li dR<0.1')
	plt.semilogy(tpr3,fpr3,'b-',label='2 Cones')
	plt.semilogy(tpr4,fpr4,'g-',label='Li dR<0.04 & L4-L1')
	plt.semilogy(tpr5,fpr5,'g:',label='Li dR<0.1 & L4/L1')
	plt.semilogy(tpr7,fpr7,'y-',label='2 Cones & L4-L1 & L4/L1')
	plt.semilogy(tpr8,fpr8,'m-',label='5 Cones')
	plt.semilogy(tpr9,fpr9,'m:',label='5 Cones Convolutional')
	plt.semilogy(tpr_csv,fpr_csv,'k-',label='CSV')
	plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	plt.title("ROC-Curves")
	plt.legend(loc=4)
	plt.savefig("ANN_data/ANNs_ROC2_log.png")
	plt.show()

def Plot_ROCs(title,tpr_no_pT, fpr_no_pT, tpr_with_pT, fpr_with_pT, tpr_with_PV, fpr_with_PV, tpr_csv, fpr_csv, x_min=0.1, x_max=0.8, y_min=0.005,folder=""):
	plt.figure("ROC_"+title)
	plt.clf()
	plt.xlim(xmin=x_min, xmax=x_max)
	plt.ylim(ymin=y_min)
	plt.semilogy(tpr_no_pT,fpr_no_pT,'r',label='no_pT')
	plt.semilogy(tpr_with_pT,fpr_with_pT,'g',label='with_pT')
	if tpr_with_PV != None and fpr_with_PV != None:
		plt.semilogy(tpr_with_PV,fpr_with_PV,'b',label='with_PV')
	plt.semilogy(tpr_csv,fpr_csv,'k-',label='CSV')
	plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	plt.title("ROC-Curves_"+title)
	plt.legend(loc=4)
	plt.savefig("Submitted_Models/"+folder+"ROCS_"+title+".png")
	print "figure saved as Submitted_Models/"+folder+"ROCS_"+title+".png"

def Plot_General_ROCs(title,datalist, x_min=0.1, x_max=0.8, y_min=0.005,folder="",DrawTitle=False):
	"""datalist should be list of tuples: (tpr, fpr, label,color(optional))"""
	c = ['r', 'g', 'b', 'brown', 'orange']
	plt.figure("ROC_"+title)
	plt.clf()
	plt.xlim(xmin=x_min, xmax=x_max)
	plt.ylim(ymin=y_min)
	for n,entry in enumerate(datalist):
		if len(entry)>3:
			plt.semilogy(entry[0],entry[1],color=entry[3],label=entry[2])
		else:
			plt.semilogy(entry[0],entry[1],color=c[n],label=entry[2])
	plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	if DrawTitle: plt.title("ROC-Curves_"+title)
	plt.legend(loc=4)
	plt.savefig("Submitted_Models/"+folder+"ROCS_"+title+".png")
	print "figure saved as Submitted_Models/"+folder+"ROCS_"+title+".png"


def Compare_GRID_ANNs():
	fpr1_csv, tpr1_csv, t_ = roc_curve(np.load('Submitted_Models/data/noPU_both/test_y.npy'),np.load('Submitted_Models/data/noPU_both/test_CSV.npy'))
	fpr3_csv, tpr3_csv, t_ = roc_curve(np.load('Submitted_Models/data/noPU_4TeV/test_y.npy'),np.load('Submitted_Models/data/noPU_4TeV/test_CSV.npy'))
	fpr5_csv, tpr5_csv, t_ = roc_curve(np.load('Submitted_Models/data/withPU_both/test_y.npy'),np.load('Submitted_Models/data/withPU_both/test_CSV.npy'))
	fpr8_csv, tpr8_csv, t_ = roc_curve(np.load('Submitted_Models/data/withPU_4TeV/test_y.npy'),np.load('Submitted_Models/data/withPU_4TeV/test_CSV.npy'))

	tpr1,fpr1 = np.loadtxt("Submitted_Models/efficiencies_noPU_both.csv",delimiter=',')
	tpr2,fpr2 = np.loadtxt("Submitted_Models/efficiencies_noPU_both_withPT.csv",delimiter=',')
	tpr3,fpr3 = np.loadtxt("Submitted_Models/efficiencies_noPU_4TeV.csv",delimiter=',')
	tpr4,fpr4 = np.loadtxt("Submitted_Models/efficiencies_noPU_4TeV_withPT.csv",delimiter=',')
	tpr5,fpr5 = np.loadtxt("Submitted_Models/efficiencies_withPU_both.csv",delimiter=',')
	tpr6,fpr6 = np.loadtxt("Submitted_Models/efficiencies_withPU_both_withPT.csv",delimiter=',')
	tpr7,fpr7 = np.loadtxt("Submitted_Models/efficiencies_withPU_both_withPV.csv",delimiter=',')
	tpr8,fpr8 = np.loadtxt("Submitted_Models/efficiencies_withPU_4TeV.csv",delimiter=',')
	tpr9,fpr9 = np.loadtxt("Submitted_Models/efficiencies_withPU_4TeV_withPT.csv",delimiter=',')
	tpr10,fpr10 = np.loadtxt("Submitted_Models/efficiencies_withPU_4TeV_withPV.csv",delimiter=',')

	Plot_ROCs("noPU_both", tpr1, fpr1, tpr2, fpr2, None, None, tpr1_csv, fpr1_csv)
	Plot_ROCs("noPU_4TeV", tpr3, fpr3, tpr4, fpr4, None, None, tpr3_csv, fpr3_csv)
	Plot_ROCs("PU_both", tpr5, fpr5, tpr6, fpr6, tpr7, fpr7, tpr5_csv, fpr5_csv)
	Plot_ROCs("PU_4TeV", tpr8, fpr8, tpr9, fpr9, tpr10, fpr10, tpr8_csv, fpr8_csv)

def Compare_GRID_ANNs_ROC1():
	fpr_csv, tpr_csv, t_ = roc_curve(np.load('Submitted_Models/data/noPU_both/test_y.npy'),np.load('Submitted_Models/data/noPU_both/test_CSV.npy'))

	tpr1,fpr1 = np.loadtxt("Submitted_Models/efficiencies_Li004.csv",delimiter=',')
	tpr2,fpr2 = np.loadtxt("Submitted_Models/efficiencies_Li01.csv",delimiter=',')
	tpr3,fpr3 = np.loadtxt("Submitted_Models/efficiencies_twoCones.csv",delimiter=',')
	tpr4,fpr4 = np.loadtxt("Submitted_Models/efficiencies_fiveCones.csv",delimiter=',')

	Plot_General_ROCs('ANNROCs1',[(tpr1,fpr1,r'Li $\Delta R<0.04$'),(tpr2,fpr2,r'Li $\Delta R<0.1$'),(tpr3,fpr3,r'two cones'),(tpr4,fpr4,r'five cones'),(tpr_csv,fpr_csv,r'CSV','k')])

def Compare_GRID_ANNs_ROC2():
	fpr_csv, tpr_csv, t_ = roc_curve(np.load('Submitted_Models/data/noPU_both/test_y.npy'),np.load('Submitted_Models/data/noPU_both/test_CSV.npy'))

	tpr1,fpr1 = np.loadtxt("Submitted_Models/efficiencies_fiveCones.csv",delimiter=',')
	tpr2,fpr2 = np.loadtxt("Submitted_Models/efficiencies_fiveConesConv.csv",delimiter=',')
	tpr3,fpr3 = np.loadtxt("Submitted_Models/efficiencies_fiveConesFunc.csv",delimiter=',')
	tpr4,fpr4 = np.loadtxt("Submitted_Models/efficiencies_noPU_both.csv",delimiter=',')

	Plot_General_ROCs('ANNROCs2',[(tpr1,fpr1,r'simple'),(tpr2,fpr2,r'conv one branch'),(tpr3,fpr3,r'conv two branches'),(tpr4,fpr4,r'conv three branches'),(tpr_csv,fpr_csv,r'CSV','k')])

def Compare_single_dataset_ANNs_on_pT_range(data_title,feature_title, pT_min, pT_max):
	if feature_title != '':
		add = "_with"+feature_title
	else:
		add = ''
	model = kr.models.load_model("Submitted_Models/model_"+data_title+add+".h5")
	test_x = np.load('Submitted_Models/data/'+data_title+'/test_x.npy')
	test_y = np.load('Submitted_Models/data/'+data_title+'/test_y.npy')
	test_CSV = np.load('Submitted_Models/data/'+data_title+'/test_CSV.npy')
	pT = np.load('Submitted_Models/data/'+data_title+'_withPT/test_feature.npy')
	if feature_title == 'PV':
		PV = np.load('Submitted_Models/data/'+data_title+'_withPV/test_feature.npy')
	
	test_x = test_x[(pT>=pT_min)*(pT<pT_max),:]	
	test_y = test_y[(pT>=pT_min)*(pT<pT_max)]	
        test_CSV= test_CSV[(pT>=pT_min)*(pT<pT_max)]
	if feature_title == 'PV':
		PV = PV[(pT>=pT_min)*(pT<pT_max)]	
	pT = pT[(pT>=pT_min)*(pT<pT_max)]
	
	if feature_title == '':
        	pred_y = model.predict(ANN_functional_shape(test_x))
	elif feature_title == 'PT':
		pred_y = model.predict(ANN_functional_shape(test_x)+[pT/200])
	elif feature_title == 'PV':
		pred_y = model.predict(ANN_functional_shape(test_x)+[PV/10])
	
	fpr, tpr, th = roc_curve(test_y,pred_y)	
	fpr_csv, tpr_csv, th_csv = roc_curve(test_y,test_CSV)
	return (fpr, tpr, fpr_csv, tpr_csv)

def Compare_GRID_ANNs_on_pT_range(pT_min,pT_max):

	folder = "pT_bins/"
	
	fpr1, tpr1, fpr1_csv, tpr1_csv = Compare_single_dataset_ANNs_on_pT_range("noPU_both",'', pT_min, pT_max)
	fpr2, tpr2, fpr1_csv, tpr1_csv = Compare_single_dataset_ANNs_on_pT_range("noPU_both",'PT', pT_min, pT_max)
	fpr3, tpr3, fpr3_csv, tpr3_csv = Compare_single_dataset_ANNs_on_pT_range("noPU_4TeV",'', pT_min, pT_max)
	fpr4, tpr4, fpr3_csv, tpr3_csv = Compare_single_dataset_ANNs_on_pT_range("noPU_4TeV",'PT', pT_min, pT_max)
	fpr5, tpr5, fpr5_csv, tpr5_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_both",'', pT_min, pT_max)
	fpr6, tpr6, fpr5_csv, tpr5_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_both",'PT', pT_min, pT_max)
	fpr7, tpr7, fpr5_csv, tpr5_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_both",'PV', pT_min, pT_max)
	fpr8, tpr8, fpr8_csv, tpr8_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_4TeV",'', pT_min, pT_max)
	fpr9, tpr9, fpr8_csv, tpr8_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_4TeV",'PT', pT_min, pT_max)
	fpr10, tpr10, fpr8_csv, tpr8_csv = Compare_single_dataset_ANNs_on_pT_range("withPU_4TeV",'PV', pT_min, pT_max)

	add = "_pT_{}_{}".format(pT_min,pT_max)

	Plot_ROCs("noPU_both"+add, tpr1, fpr1, tpr2, fpr2, None, None, tpr1_csv, fpr1_csv,folder=folder)
	#Plot_ROCs("noPU_4TeV"+add, tpr3, fpr3, tpr4, fpr4, None, None, tpr3_csv, fpr3_csv,folder=folder)
	Plot_ROCs("PU_both"+add, tpr5, fpr5, tpr6, fpr6, tpr7, fpr7, tpr5_csv, fpr5_csv,folder=folder)
	#Plot_ROCs("PU_4TeV"+add, tpr8, fpr8, tpr9, fpr9, tpr10, fpr10, tpr8_csv, fpr8_csv,folder=folder)


if __name__ == "__main__":

	#Compare_GRID_ANNs()
	#pT_bins = [0,1200,1800,2500]
	
	'''
	pT_bins = [0,1200,2500]

	for n in range(len(pT_bins)-1):
		Compare_GRID_ANNs_on_pT_range(pT_bins[n],pT_bins[n+1])
	'''
	Compare_GRID_ANNs_ROC1()
	Compare_GRID_ANNs_ROC2()	
	
	#initialize model
	
	#model = kr.models.load_model("ANN_model_complex_PU_pT2_4TeVonly.h5")
	#model = build_complex_functional_model(Additional_Input=True)
	#model = build_conv_model()
	#print "model initialized"
	
	
	#preprocess data
	
	#New_Data('PU2_4TeVonly')
	#New_Data_4_Submit()
	
	'''
	test_x  =np.load("ANN_data/test_x_PU2_4TeVonly.npy")
	test_y  =np.load("ANN_data/test_y_PU2_4TeVonly.npy")
	test_CSV=np.load("ANN_data/test_CSV_PU2_4TeVonly.npy")
	test_pT =np.load("ANN_data/test_pT_PU2_4TeVonly.npy")
	train_x =np.load("ANN_data/train_x_PU2_4TeVonly.npy")
	train_y =np.load("ANN_data/train_y_PU2_4TeVonly.npy")
	train_pT =np.load("ANN_data/train_pT_PU2_4TeVonly.npy")
	'''
	
	'''
	#Convolutional ANN
	m_train = train_x.shape[0]
	m_test = test_x.shape[0]
	test_x = np.reshape(test_x.flatten(),(-1,5,4,1))
	train_x = np.reshape(train_x.flatten(),(-1,5,4,1))
	'''
	'''
	#functional ANN
	m_train = train_x.shape[0]
	m_test = test_x.shape[0]
	
	test_x_Li_Lj=discriminants(test_x)
	train_x_Li_Lj=discriminants(train_x)
	test_x_Li = np.reshape(test_x[:,:20].flatten(),(-1,5,4,1))
	train_x_Li = np.reshape(train_x[:,:20].flatten(),(-1,5,4,1))
	test_x_PU = test_x[:,-1]
	train_x_PU = train_x[:,-1]
	
	print "data loaded"
	'''
	'''
	#sample_weights for PU
	from WeightFunctions import reweight
	weights = reweight(train_x_PU, train_y, "PU")
	
	test_x_PU = test_x_PU/10
	train_x_PU = train_x_PU/10
	'''
	'''
	#sample_weights for jet_pT
	from WeightFunctions import reweight
	weights = reweight(train_pT, train_y, "jet_pT")
	
	train_pT = train_pT/200
	test_pT = test_pT/200
	
	#Compare_ANNs()
	
	#history = model.fit(train_x_Li, train_y, epochs=90, batch_size=32, validation_data=(test_x_Li,test_y))
	#history = model.fit([train_x_Li,train_x_Li_Lj], train_y, epochs=60, batch_size=32, validation_data=([test_x_Li,test_x_Li_Lj],test_y))
	#history = model.fit([train_x_Li,train_x_Li_Lj,train_x_PU], train_y, epochs=90, batch_size=32, validation_data=([test_x_Li,test_x_Li_Lj,test_x_PU],test_y),sample_weight=weights)
	history = model.fit([train_x_Li,train_x_Li_Lj,train_pT], train_y, epochs=20, batch_size=32, validation_data=([test_x_Li,test_x_Li_Lj,test_pT],test_y),sample_weight=weights)
	
	model.save("ANN_model_complex_PU_pT2_4TeVonly.h5")
	
	plt.figure('History')
	plt.plot(history.history['loss'],label='training')
	plt.plot(history.history['val_loss'],label='test')
	plt.xlabel('epoch')
	plt.ylabel('loss')
	plt.legend(loc=2)
	
	
	pred_y = model.predict([test_x_Li,test_x_Li_Lj,test_pT])
	fpr, tpr, thresholds = roc_curve(test_y,pred_y)
	fpr_csv, tpr_csv, thresholds_csv = roc_curve(test_y,test_CSV)
	
	plt.figure('ROC')
	plt.semilogy(tpr,fpr,'r-',label='model')
	plt.semilogy(tpr_csv,fpr_csv,'b-',label='CSV')
	plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.legend(loc=3)
	plt.show()
	
	
	csv_file = open("ANN_data/complex_PU_pT2_4TeVonly.csv","wb")
	writer = csv.writer(csv_file)
	writer.writerow(tpr)
	writer.writerow(fpr)
	csv_file.close()
	'''
	'''
	tpr_PU,fpr_PU = np.loadtxt("ANN_data/complex_PU.csv",delimiter=',')
	fpr_csv_PU, tpr_csv_PU, thresholds_csv_PU = roc_curve(test_y,test_CSV)
	
	tpr_PU2,fpr_PU2 = np.loadtxt("ANN_data/complex_PU_new.csv",delimiter=',')
	
	tpr9,fpr9 = np.loadtxt("ANN_data/5Cone_conv.csv",delimiter=',')
	test_CSV2=np.load("ANN_data/test_CSV.npy")
	test_y2  =np.load("ANN_data/test_y.npy")
	fpr_csv, tpr_csv, thresholds_csv = roc_curve(test_y2,test_CSV2)
	
	plt.figure("ROC")
	#plt.semilogy(tpr,fpr,label='complex_PU')
	plt.semilogy(tpr9,fpr9,'r-',label='best_noPU')
	plt.semilogy(tpr_csv,fpr_csv,'b-',label='CSV_noPU')
	plt.semilogy(tpr_PU,fpr_PU,'r-.',label='best_PU')
	plt.semilogy(tpr_PU2,fpr_PU2,'g-.',label='best_PU_using_PV')
	plt.semilogy(tpr_csv_PU,fpr_csv_PU,'b-.',label='CSV_PU')
	plt.semilogy([0,1],[0.1,0.1],'k:',label="10% mistag")
	plt.xlabel(r"$\epsilon$_signal")
	plt.ylabel(r"$\epsilon$_background")
	plt.title("ROC-Curves")
	plt.legend(loc=3)
	plt.savefig("ROC/ANN_ROC_PU.png")
	
	plt.show()
	'''








