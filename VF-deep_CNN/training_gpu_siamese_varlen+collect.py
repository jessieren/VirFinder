#source /usr/usc/cuda/default/setup.sh
#source /usr/usc/cuDNN/default/setup.sh

#export LD_LIBRARY_PATH=/usr/usc/cuDNN/default//lib64:$LD_LIBRARY_PATH
#export CPATH=/usr/usc/cuDNN/default/include:$CPATH
#export LIBRARY_PATH=/usr/usc/cuDNN/default/lib64:$LIBRARY_PATH

###################################################################

def get_output(input_layer, hidden_layers):
    output = input_layer
    for hidden_layer in hidden_layers:
        output = hidden_layer(output)
    return output


#import multiprocessing
import os, sys, glob, path
#os.environ['THEANO_FLAGS'] = "floatX=float32,openmp=True" 
#os.environ['THEANO_FLAGS'] = "mode=FAST_RUN,device=gpu0,floatX=float32" 
#os.environ['OMP_NUM_THREADS'] = str(multiprocessing.cpu_count())


import optparse
prog_base = os.path.split(sys.argv[0])[1]

parser = optparse.OptionParser()
parser.add_option("-l", "--len", action = "store", type = int, dest = "contigLength",
									help = "contig Length")
parser.add_option("-i", "--in", action = "store", type = "string", dest = "inDir",
									default='./', help = "input directory")
parser.add_option("-o", "--out", action = "store", type = "string", dest = "outDir",
									default='./', help = "output directory")
parser.add_option("-f", "--fLen1", action = "store", type = int, dest = "filter_len1",
									help = " filter length layer 1 ")
parser.add_option("-n", "--fNum1", action = "store", type = int, dest = "nb_filter1",
									default=0, help = "number of filter layer 1")
parser.add_option("-d", "--dense", action = "store", type = int, dest = "nb_dense",
									default=0, help = "number of dense states")
parser.add_option("-e", "--epochs", action = "store", type = int, dest = "epochs",
									default=0, help = "number of epochs")
#parser.add_option("-u", "--continue", action = "store", type = "string", dest = "onlyComputeMeasure",
#									default=0, help = "kmer count is ready only compute measures? 1 for yes, 0 for no")
#parser.add_option("-k", "--kLen", action = "store", type = "string", dest = "kLen",
#									help = "the length of k-tuple")

(options, args) = parser.parse_args()
if (options.contigLength is None or
		options.filter_len1 is None or
    options.nb_filter1 is None or options.nb_dense is None ) :
	sys.stderr.write(prog_base + ": ERROR: missing required command-line argument")
	filelog.write(prog_base + ": ERROR: missing required command-line argument")
	parser.print_help()
	sys.exit(0)


contigLength = options.contigLength
filter_len1 = options.filter_len1
nb_filter1 = options.nb_filter1
nb_dense = options.nb_dense
inDir = options.inDir
outDir = options.outDir
epochs = options.epochs
#print('contigLength='+str(contigLength))

#import os
#os.environ["THEANO_FLAGS"] = "mode=FAST_RUN,device=gpu0,floatX=float32"



import numpy as np
import keras
from keras.models import Sequential, load_model, Model
from keras.layers import Dense, Dropout, Activation, Flatten, Layer, Input, Conv1D, MaxPooling1D, GlobalMaxPooling1D
from keras.layers.merge import Average
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import Adagrad, Adam
import h5py
channel_num = 4



#contigLength = 1000
#contigLengthk = str(int(contigLength/1000))
contigLengthk = contigLength/1000
if contigLengthk.is_integer() :
    contigLengthk = int(contigLengthk)
contigLengthk = str(contigLengthk)

#outDir = '/staging/fs3/renj/'
print("...loading data...")

# Y_1k_shuf_val128186_code.npy
Xname_tr_fw = [ x for x in os.listdir(inDir) if 'X_'+str(contigLengthk)+'k_shuf_trfw' in x and '_code.npy' in x ][0]
Xname_tr_bw = Xname_tr_fw.replace('_trfw', '_trbw')
Yname_tr = Xname_tr_fw.replace('_trfw', '_tr').replace('X_', 'Y_')
print(Yname_tr)
Y_tr_shuf = np.load( os.path.join(inDir, Yname_tr ) )
X_tr_shuf_fw = np.load( os.path.join(inDir, Xname_tr_fw ) )
X_tr_shuf_bw = np.load( os.path.join(inDir, Xname_tr_bw ) )

Xname_val_fw = [ x for x in os.listdir(inDir) if 'X_'+str(contigLengthk)+'k_shuf_valfw' in x and '_code.npy' in x ][0]
Xname_val_bw = Xname_val_fw.replace('_valfw', '_valbw')
Yname_val = Xname_val_fw.replace('_valfw', '_val').replace('X_', 'Y_')
print(Yname_val)
Y_val_shuf = np.load( os.path.join(inDir,  Yname_val ) )
X_val_shuf_fw = np.load( os.path.join(inDir, Xname_val_fw ) )
X_val_shuf_bw = np.load( os.path.join(inDir, Xname_val_bw ) )




# parameters
#filter_len1 = 6
#filter_len2 = 10
#nb_filter1 = 1000
#nb_filter2 = 100
#nb_hidden = 100
POOL_FACTOR = 1
dropout_cnn = 0.1
dropout_pool = 0.1
dropout_dense = 0.1
learningrate = 0.001
#LR = 0.01/2
batch_size=int(X_tr_shuf_fw.shape[0]/(2000*1000/contigLength)) ## smaller batch size can reduce memory
#batch_size=10 ## smaller batch size can reduce memory
#epochs=1

pool_len1 = int((contigLength-filter_len1+1)/POOL_FACTOR)
#pool_len2 = int((seq_len-filter_len2+1)/POOL_FACTOR)
#pool_len1 = 10


modPattern = 'model_siamese_varlen_'+contigLengthk+'k_fl'+str(filter_len1)+'_fn'+str(nb_filter1)+'_dn'+str(nb_dense)
modName = os.path.join( outDir, modPattern +'_ep'+str(epochs)+'.h5')
checkpointer = keras.callbacks.ModelCheckpoint(filepath=modName, verbose=1, save_best_only=True)
earlystopper = keras.callbacks.EarlyStopping(monitor='val_acc', min_delta=0.0001, patience=5, verbose=1)


##### build model #####
print("...building model...")
## if model exists
if os.path.isfile(modName):
  model = load_model(modName)
  print("...model exists...")
else :
  ## siamese
  forward_input = Input(shape=(None, channel_num))
  reverse_input = Input(shape=(None, channel_num))
  hidden_layers = [
      Conv1D(filters = nb_filter1, kernel_size = filter_len1, activation='relu'),
      #Dropout(dropout_cnn),
      GlobalMaxPooling1D(),
      # https://github.com/fchollet/keras/issues/1920
      # https://stats.stackexchange.com/questions/257321/what-is-global-max-pooling-layer-and-what-is-its-advantage-over-maxpooling-layer
      #MaxPooling1D(pool_len1),
      Dropout(dropout_pool),
      #Flatten(),
      Dense(nb_dense, activation='relu'),
      Dropout(dropout_dense),
      Dense(1, activation='sigmoid')
  ]
  forward_output = get_output(forward_input, hidden_layers)     
  reverse_output = get_output(reverse_input, hidden_layers)
  #output = merge([forward_output, reverse_output], mode='ave')
  output = Average()([forward_output, reverse_output])
  model = Model(inputs=[forward_input, reverse_input], outputs=output)
  model.compile(Adam(lr=learningrate), 'binary_crossentropy', metrics=['accuracy'])
  #model.summary()


print("...fitting model...")
print(contigLengthk+'k_fl'+str(filter_len1)+'_fn'+str(nb_filter1)+'_dn'+str(nb_dense)+'_ep'+str(epochs))
model.fit(x = [X_tr_shuf_fw, X_tr_shuf_bw], y = Y_tr_shuf, \
            batch_size=batch_size, epochs=epochs, verbose=1, \
            validation_data=([X_val_shuf_fw, X_val_shuf_bw], Y_val_shuf), \
            callbacks=[checkpointer, earlystopper])
            
            
            
### produce AUC ###
            
import sklearn
from sklearn.metrics import roc_auc_score 


## train data
type = 'tr'
print("...predicting "+type+"...\n")

Y_shuf_pred = model.predict([X_tr_shuf_fw, X_tr_shuf_bw], batch_size=1)
auc = sklearn.metrics.roc_auc_score(Y_tr_shuf, Y_shuf_pred)
print('auc_'+type+'='+str(auc)+'\n')
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_pred.txt'), np.transpose(Y_shuf_pred))
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_true.txt'), np.transpose(Y_tr_shuf))

del Y_tr_shuf
del X_tr_shuf_fw
del X_tr_shuf_bw
del Y_shuf_pred

## val data
type = 'val'
print("...predicting "+type+"...\n")

Y_shuf_pred = model.predict([X_val_shuf_fw, X_val_shuf_bw], batch_size=1)
auc = sklearn.metrics.roc_auc_score(Y_val_shuf, Y_shuf_pred)
print('auc_'+type+'='+str(auc)+'\n')
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_pred.txt'), np.transpose(Y_shuf_pred))
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_true.txt'), np.transpose(Y_val_shuf))

del Y_val_shuf
del X_val_shuf_fw
del X_val_shuf_bw
del Y_shuf_pred


## test data
type = 'te'
print("...loading "+type+" data...")

Xname_fw = [ x for x in os.listdir(inDir) if 'X_'+str(contigLengthk)+'k_shuf_'+type+'fw' in x and '_code.npy' in x ][0]
Xname_bw = Xname_fw.replace('_'+type+'fw', '_'+type+'bw')
Yname = Xname_fw.replace('_'+type+'fw', '_'+type ).replace('X_', 'Y_')
print(Yname)
Y_shuf = np.load( os.path.join(inDir, Yname ) )
X_shuf_fw = np.load( os.path.join(inDir, Xname_fw ) )
X_shuf_bw = np.load( os.path.join(inDir, Xname_bw ) )

print("...predicting "+type+"...\n")

Y_shuf_pred = model.predict([X_shuf_fw, X_shuf_bw], batch_size=1)
auc = sklearn.metrics.roc_auc_score(Y_shuf, Y_shuf_pred)
print('auc_'+type+'='+str(auc)+'\n')
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_pred.txt'), np.transpose(Y_shuf_pred))
np.savetxt(os.path.join(outDir, modPattern + '_' + type + 'fw_Y_true.txt'), np.transpose(Y_shuf))

del Y_shuf
del X_shuf_fw
del X_shuf_bw



