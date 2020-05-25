#!/usr/bin/env python
"""
Script for predict single cell TF binding using pre-trained model.

Use `scFAN_predict.py -h` to see an auto-generated description of advanced options.
python scFAN_predict.py -i /data2/fly/PBMCs/raw_data/gz_files/new_folder/LMPP -scindir /data2/fly/scFAN_data/new_folder_cisTopics_PBMC_agg -moname LMPP -pb True -oc multiTask_H1hESC_add_ATAC_moreTFs multiTask_GM12878_add_ATAC_moreTFs multiTask_K562_ATAC_more_chipdata
"""
import utils_scFAN as utils
import numpy as np

# Standard library imports
import sys
import os
import errno
import argparse
import pickle
import pdb
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import precision_recall_curve
from collections import Counter
import ast
import commands
#import tensorflow as tf
#config = tf.ConfigProto(device_count={'gpu':3})
#config.gpu_options.allow_growth=True
#session = tf.Session(config=config)

def test(model_dir,datagen_test):
    model_tfs_list = []
    model_predicts_list = []
    for model_dir_item in model_dir:
        model_tfs, model_bigwig_names, features, model = utils.load_model(model_dir_item)
        #pdb.set_trace()
        model_predicts = model.predict_generator(datagen_test, val_samples=1+len(datagen_test)/100, pickle_safe=True,verbose=1)  ## old_version
    #model_predicts = model.predict_generator(datagen_test, steps=1+len(datagen_test)/100, use_multiprocessing=True,verbose=1)
        model_tfs_list.append(model_tfs)
        model_predicts_list.append(model_predicts)
    return model_predicts_list,model_tfs_list

def train(datagen_train, datagen_valid, model, epochs, patience, learningrate, output_dir):
    from keras.callbacks import ModelCheckpoint, EarlyStopping
    from keras.optimizers import Adam

    print 'Compiling model'
    model.compile(Adam(lr=learningrate), 'binary_crossentropy', metrics=['accuracy'])

    model.summary()

    print 'Running at most', str(epochs), 'epochs'

    checkpointer = ModelCheckpoint(filepath=output_dir + '/best_model.hdf5',
                                   verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=patience, verbose=1)

    train_samples_per_epoch = len(datagen_train)/epochs/utils.batch_size*utils.batch_size
    history = model.fit_generator(datagen_train, samples_per_epoch=train_samples_per_epoch,
                                  nb_epoch=epochs, validation_data=datagen_valid,
                                  nb_val_samples=len(datagen_valid),
                                  callbacks=[checkpointer, earlystopper],
                                  pickle_safe=True)

    print 'Saving final model'
    model.save_weights(output_dir + '/final_model.hdf5', overwrite=True)

    print 'Saving history'
    history_file = open(output_dir + '/history.pkl', 'wb')
    pickle.dump(history.history, history_file)
    history_file.close()


def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Train model.",
        epilog='\n'.join(__doc__.strip().split('\n')[1:]).strip(),
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--inputdirs', '-i', type=str, required=True, nargs='+',
                        help='Folders containing data.')    
    parser.add_argument('--input_scATAC_dir', '-scindir', type=str, required=True, nargs='+',
                        help='Folders containing scATAC-seq data.')    
    parser.add_argument('--epochs', '-e', type=int, required=False,
                        default=100,
                        help='Epochs to train (default: 100).')
    parser.add_argument('--patience', '-ep', type=int, required=False,
                        default=20,
                        help='Number of epochs with no improvement after which training will be stopped (default: 20).')
    parser.add_argument('--learningrate', '-lr', type=float, required=False,
                        default=0.001,
                        help='Learning rate (default: 0.001).')
    parser.add_argument('--negatives', '-n', type=int, required=False,
                        default=1,
                        help='Number of negative samples per each positive sample (default: 1).')
    parser.add_argument('--seqlen', '-L', type=int, required=False,
                        default=1000,
                        help='Length of sequence input (default: 1000).')
    parser.add_argument('--bigwig_weight', '-bgw', type=int, required=False,
                        default=0.75,
                        help='bigwig file weights.')
    parser.add_argument('--dense', '-d', type=int, required=False,
                        default=128,
                        help='Number of dense units in model (default: 128).')
    parser.add_argument('--dropout', '-p', type=float, required=False,
                        default=0.5,
                        help='Dropout rate between the LSTM and dense layers (default: 0.5).')
    parser.add_argument('--seed', '-s', type=int, required=False,
                        default=420,
                        help='Random seed for consistency (default: 420).')
    parser.add_argument('--factor', '-f', type=str, required=False,
                        default=None,
                        help='The transcription factor to train. If not specified, multi-task training is used instead.')
    parser.add_argument('--meta', '-m', action='store_true',
                        help='Meta flag. If used, model will use metadata features.')
    parser.add_argument('--gencode', '-g', action='store_true',
                        help='GENCODE flag. If used, model will incorporate CpG island and gene annotation features.')
    parser.add_argument('--process_batch', '-pb', type=ast.literal_eval,
                        help='whether process situation where batch exists, default is true')
    parser.add_argument('--motif', '-mo', action='store_true',
                        help='Motif flag. If used, will inject canonical motif and its RC as model model weights (if available).')
    parser.add_argument('--motif_name', '-moname', type=str,required=True,default='GM',
                        help='Motif flag. If used, will inject canonical motif and its RC as model model weights (if available).')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-o', '--outputdir', type=str,
                       help='The output directory. Causes error if the directory already exists.')
    group.add_argument('-oc', '--outputdirc', type=str,nargs='+',
                       help='The output directory. Will overwrite if directory already exists.')
    return parser


def main():
    """
    The main executable function
    """
    parser = make_argument_parser()
    args = parser.parse_args()

    input_dirs = args.inputdirs
    input_scATAC_dir = args.input_scATAC_dir
    tf = args.factor
    epochs = args.epochs
    patience = args.patience
    process_batch = args.process_batch
    learningrate = args.learningrate
    seed = args.seed
    utils.set_seed(seed)
    dropout_rate = args.dropout
    L = args.seqlen
    utils.L = L
    negatives = args.negatives
    assert negatives > 0
    meta = args.meta
    gencode = args.gencode
    motif = args.motif
    motif_name = args.motif_name

    num_dense = args.dense
 
    features = ['bigwig']    
    #pdb.set_trace()


    if args.outputdir is None:
        clobber = True
        output_dir = args.outputdirc
    else:
        clobber = False
        output_dir = args.outputdir

    try:  # adapted from dreme.py by T. Bailey
        os.makedirs(output_dir[0])
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not clobber:
                print >> sys.stderr, ('pretrained model directory (%s) already exists '
                                      'but you specified not to clobber it') % output_dir[0]
                sys.exit(1)
            else:
                print >> sys.stderr, ('pretrained model directory (%s) already exists '
                                      'so it will be clobbered') % output_dir[0]

    print 'Loading genome'
    genome = utils.load_genome()
    print 'Loading ChIP labels'
    assert len(input_dirs) == 1 # multi-task training only supports one cell line
    input_dir = input_dirs[0]
    tfs, positive_windows_list = \
            utils.load_chip_multiTask_multiple(input_dir,process_batch)
    num_tfs = len(tfs)
    print 'Loading bigWig data'
    bigwig_names, bigwig_files_list = utils.load_bigwigs(input_dirs)
    num_bigwigs = len(bigwig_names)
    chip_name_file = np.loadtxt(input_dir + '/chip.txt',dtype=str)
    big_wig_list = [item.split('_')[0]+'.bw' for item in chip_name_file[:,0]]
    model_predicts_list = []
    test_predict_all_TFs_H1 = []
    bb = []
    test_predict_all_TFs_GM = []
    test_predict_all_TFs_K562 = []
    TF_score_H1 = []
    TF_score_GM = []
    TF_score_K562 = []
    for num_pos,positive_windows in enumerate(positive_windows_list):
        _, datagen_bed = utils.load_bed_data_sc(genome, positive_windows, False, False, input_dir, False, big_wig_list,num_pos,input_scATAC_dir, chrom=None)
        #pdb.set_trace()
        
        print "%d sample...in %d "%(num_pos+1,len(positive_windows_list))
        model_predicts_list,model_tfs_list = test(output_dir,datagen_bed)
        try:
            os.stat('%s/scFAN_predict_using_H1'%(input_dirs[0]))
        except:
            os.mkdir('%s/scFAN_predict_using_H1'%(input_dirs[0]))
            os.mkdir('%s/scFAN_predict_using_GM'%(input_dirs[0]))
            os.mkdir('%s/scFAN_predict_using_K562'%(input_dirs[0]))
        np.save('%s/scFAN_predict_using_H1/%s_data'%(input_dirs[0],str(num_pos)),model_predicts_list[0])
        np.save('%s/scFAN_predict_using_GM/%s_data'%(input_dirs[0],str(num_pos)),model_predicts_list[1])
        np.save('%s/scFAN_predict_using_K562/%s_data'%(input_dirs[0],str(num_pos)),model_predicts_list[2])
        print 'calculating TF activity score for cell %s'%(str(num_pos))
        for index in range(model_predicts_list[0].shape[0]):
            test_predict_all_TFs_H1 = test_predict_all_TFs_H1 + [model_tfs_list[0][item] for item in np.argsort(model_predicts_list[0][index])[::-1][:2]]
            test_predict_all_TFs_GM = test_predict_all_TFs_GM + [model_tfs_list[1][item] for item in np.argsort(model_predicts_list[1][index])[::-1][:2]]
            test_predict_all_TFs_K562 = test_predict_all_TFs_K562 + [model_tfs_list[2][item] for item in np.argsort(model_predicts_list[2][index])[::-1][:2]]
        TF_score_H1.append([Counter(test_predict_all_TFs_H1)[item]/float(sum(Counter(test_predict_all_TFs_H1).values())) for item in model_tfs_list[0]])
        TF_score_GM.append([Counter(test_predict_all_TFs_GM)[item]/float(sum(Counter(test_predict_all_TFs_GM).values())) for item in model_tfs_list[1]])
        TF_score_K562.append([Counter(test_predict_all_TFs_K562)[item]/float(sum(Counter(test_predict_all_TFs_K562).values())) for item in model_tfs_list[2]])
        test_predict_all_TFs_H1 = []
        test_predict_all_TFs_GM = []
        test_predict_all_TFs_K562 = []
        
    #pdb.set_trace()
    TF_score_H1 = np.array(TF_score_H1)
    TF_score_GM = np.array(TF_score_GM)
    TF_score_K562 = np.array(TF_score_K562)
    print 'Done calculating TF activity score...saving results!!!'
    np.save('%s/TF_activity_score_%s_pretrained_model_H1'%(input_dirs[0],motif_name),TF_score_H1)
    np.save('%s/TF_activity_score_%s_pretrained_model_GM'%(input_dirs[0],motif_name),TF_score_GM)
    np.save('%s/TF_activity_score_%s_pretrained_model_K562'%(input_dirs[0],motif_name),TF_score_K562)
    '''
    print 'Cleaning tmp files...'
    t1 = commands.getoutput('rm -f %s/scFAN_predict_using_H1/*.npy'%(input_dirs[0]))
    t2 = commands.getoutput('rm -f %s/scFAN_predict_using_GM/*.npy'%(input_dirs[0]))
    t3 = commands.getoutput('rm -f %s/scFAN_predict_using_K562/*.npy'%(input_dirs[0]))
    '''
if __name__ == '__main__':
    """
    See module-level docstring for a description of the script.
    """
    main()
