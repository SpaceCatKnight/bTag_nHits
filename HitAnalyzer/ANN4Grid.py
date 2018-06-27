import ROOT as rt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import keras as kr
import pickle
import csv
from sklearn.metrics import roc_curve
import sys

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

if __name__ == "__main__":

	if len(sys.argv) > 1: title = sys.argv[1]
        if len(sys.argv) > 2: data_path = sys.argv[2]
	if len(sys.argv) > 3: addFeature = sys.argv[3]
	if len(sys.argv) > 4: nEpochs = int(sys.argv[4])
	if len(sys.argv) > 5: modenr = int(sys.argv[5])

	#initialize model
	if addFeature == "No":
		model = build_complex_functional_model(Additional_Input=False)
	else:
		model = build_complex_functional_model(Additional_Input=True)


	#preprocess data

	test_x  =np.load(data_path+"/test_x.npy")
	test_y  =np.load(data_path+"/test_y.npy")
	test_CSV=np.load(data_path+"/test_CSV.npy")
	
	train_x =np.load(data_path+"/train_x.npy")
	train_y =np.load(data_path+"/train_y.npy")

	if addFeature != "No":	
		test_feature =np.load(data_path+"/test_feature.npy")
		train_feature =np.load(data_path+"/train_feature.npy")

	m_train = train_x.shape[0]
	m_test = test_x.shape[0]
	
	test_x_Li_Lj=discriminants(test_x)
	train_x_Li_Lj=discriminants(train_x)
	test_x_Li = np.reshape(test_x[:,:20].flatten(),(-1,5,4,1))
	train_x_Li = np.reshape(train_x[:,:20].flatten(),(-1,5,4,1))
	test_x_PU = test_x[:,-1]
	train_x_PU = train_x[:,-1]

	if addFeature == 'PV':
		from WeightFunctions import reweight
		weights = reweight(train_feature, train_y, "PU", modenr)

		test_feature = test_feature/10
		train_feature = train_feature/10
	
	elif addFeature == 'pT':
		from WeightFunctions import reweight
		weights = reweight(train_feature, train_y, "jet_pT", modenr)

		train_feature = train_feature/200
		test_feature = test_feature/200


	if addFeature == "No":
		history = model.fit([train_x_Li,train_x_Li_Lj], train_y, epochs=nEpochs, verbose=2, batch_size=32, validation_data=([test_x_Li,test_x_Li_Lj],test_y))
	else:
		history = model.fit([train_x_Li,train_x_Li_Lj,train_feature], train_y, epochs=nEpochs, verbose=2, batch_size=32, validation_data=([test_x_Li,test_x_Li_Lj,test_feature],test_y),sample_weight=weights)

	model.save("Submitted_Models/model_{}.h5".format(title))

	plt.figure('History')
	plt.plot(history.history['loss'],label='training')
	plt.plot(history.history['val_loss'],label='test')
	plt.xlabel('epoch')
	plt.ylabel('loss')
	plt.legend(loc=2)
	plt.savefig("Submitted_Models/history_{}.png".format(title))

	if addFeature == "No":
		pred_y = model.predict([test_x_Li,test_x_Li_Lj])
	else:
		pred_y = model.predict([test_x_Li,test_x_Li_Lj,test_feature])
	fpr, tpr, thresholds = roc_curve(test_y,pred_y)
	fpr_csv, tpr_csv, thresholds_csv = roc_curve(test_y,test_CSV)

	plt.figure('ROC')
	plt.semilogy(tpr,fpr,'r-',label='model')
	plt.semilogy(tpr_csv,fpr_csv,'b-',label='CSV')
	plt.semilogy([0,1],[0.1,0.1],'k:',label='10% mistag')
	plt.legend(loc=3)
	plt.savefig("Submitted_Models/ROC_{}.png".format(title))

	csv_file = open("Submitted_Models/efficiencies_{}.csv".format(title),"wb")
	writer = csv.writer(csv_file)
	writer.writerow(tpr)
	writer.writerow(fpr)
	csv_file.close()






