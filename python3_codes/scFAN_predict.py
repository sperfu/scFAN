#!/usr/bin/env python
"""
Script for training model.
-k 128 -r 64 -d 256 -e 5
Use `train.py -h` to see an auto-generated description of advanced options.
python scFAN_predict.py -i data/H1-hESC/scATAC_data/gz_files_SU353-LSC -moname SU353-LSC -oc multiTask_H1hESC_add_ATAC_moreTFs
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


def test(model_dir,datagen_test):
    model_tfs, model_bigwig_names, features, model = utils.load_model(model_dir)
    #pdb.set_trace()
    model_predicts = model.predict_generator(datagen_test, val_samples=1+len(datagen_test)/100, pickle_safe=True,verbose=1)  ## old_version
    #model_predicts = model.predict_generator(datagen_test, steps=1+len(datagen_test)/100, use_multiprocessing=True,verbose=1)
    return model_predicts,model_tfs

def train(datagen_train, datagen_valid, model, epochs, patience, learningrate, output_dir):
    from keras.callbacks import ModelCheckpoint, EarlyStopping
    from keras.optimizers import Adam

    print('Compiling model')
    model.compile(Adam(lr=learningrate), 'binary_crossentropy', metrics=['accuracy'])

    model.summary()

    print('Running at most', str(epochs), 'epochs')

    checkpointer = ModelCheckpoint(filepath=output_dir + '/best_model.hdf5',
                                   verbose=1, save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=patience, verbose=1)

    train_samples_per_epoch = len(datagen_train)/epochs/utils.batch_size*utils.batch_size
    history = model.fit_generator(datagen_train, samples_per_epoch=train_samples_per_epoch,
                                  nb_epoch=epochs, validation_data=datagen_valid,
                                  nb_val_samples=len(datagen_valid),
                                  callbacks=[checkpointer, earlystopper],
                                  pickle_safe=True)

    print('Saving final model')
    model.save_weights(output_dir + '/final_model.hdf5', overwrite=True)

    print('Saving history')
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
    parser.add_argument('--validinputdirs', '-vi', type=str, required=False, nargs='+',
                        default=None,
                        help='Folder(s) of validation cell type data (optional).')
    parser.add_argument('--validchroms', '-v', type=str, required=False, nargs='+',
                        default=['chr11'],
                        help='Chromosome(s) to set aside for validation (default: chr11).')
    parser.add_argument('--testchroms', '-t', type=str, required=False, nargs='+',
                        default=['chr1', 'chr8', 'chr21'],
                        help='Chromosome(s) to set aside for testing (default: chr1, chr8, chr21).')
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
    parser.add_argument('--motifwidth', '-w', type=int, required=False,
                        default=26,
                        help='Width of the convolutional kernels (default: 26).')
    parser.add_argument('--bigwig_weight', '-bgw', type=int, required=False,
                        default=0.75,
                        help='bigwig file weights.')
    parser.add_argument('--kernels', '-k', type=int, required=False,
                        default=32,
                        help='Number of kernels in model (default: 32).')
    parser.add_argument('--recurrent', '-r', type=int, required=False,
                        default=0,
                        help='Number of LSTM units in model (default: 0). If set to 0, the BLSTM layer is simply removed. If negative, the BLSTM layer is replaced with a global max-pooling layer.')
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
    parser.add_argument('--motif', '-mo', action='store_true',
                        help='Motif flag. If used, will inject canonical motif and its RC as model model weights (if available).')
    parser.add_argument('--motif_name', '-moname', type=str,required=True,default='GM',
                        help='Motif flag. If used, will inject canonical motif and its RC as model model weights (if available).')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-o', '--outputdir', type=str,
                       help='The output directory. Causes error if the directory already exists.')
    group.add_argument('-oc', '--outputdirc', type=str,
                       help='The output directory. Will overwrite if directory already exists.')
    return parser


def main():
    """
    The main executable function
    """
    parser = make_argument_parser()
    args = parser.parse_args()

    input_dirs = args.inputdirs
    tf = args.factor
    valid_chroms = args.validchroms
    valid_input_dirs = args.validinputdirs
    test_chroms = args.testchroms
    epochs = args.epochs
    patience = args.patience
    learningrate = args.learningrate
    seed = args.seed
    utils.set_seed(seed)
    dropout_rate = args.dropout
    L = args.seqlen
    w = args.motifwidth
    utils.L = L
    utils.w = w
    utils.w2 = w/2
    negatives = args.negatives
    assert negatives > 0
    meta = args.meta
    gencode = args.gencode
    motif = args.motif
    motif_name = args.motif_name

    num_motifs = args.kernels
    num_recurrent = args.recurrent
    num_dense = args.dense
 
    features = ['bigwig']    
    #pdb.set_trace()

    if tf:
        print('Single-task training:', tf)
        
    else:
        print('Multi-task training')
        singleTask = False
        #Cannot use any metadata features
        assert not meta
        assert not gencode

    if args.outputdir is None:
        clobber = True
        output_dir = args.outputdirc
    else:
        clobber = False
        output_dir = args.outputdir

    try:  # adapted from dreme.py and train.py by T. Bailey & Daniel Quang
        os.makedirs(output_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not clobber:
                print >> sys.stderr, ('output directory (%s) already exists '
                                      'but you specified not to clobber it') % output_dir
                sys.exit(1)
            else:
                print >> sys.stderr, ('output directory (%s) already exists '
                                      'so it will be clobbered') % output_dir

    print('Loading genome')
    genome = utils.load_genome()
    if valid_input_dirs:
        print('You specified at least one validation input directory')
        assert singleTask # This option only works for single-task training
    print('Loading ChIP labels')
    if singleTask:
        num_tfs = 1
    else:
        assert len(input_dirs) == 1 # multi-task training only supports one cell line
        input_dir = input_dirs[0]
        #tfs, positive_windows, y_positive, nonnegative_regions_bed = \
        tfs, positive_windows_list, y_positive, nonnegative_regions_bed = \
            utils.load_chip_multiTask_multiple(input_dir)
            #utils.load_chip_multiTask(input_dir)
        num_tfs = len(tfs)
    print('Loading bigWig data')
    bigwig_names, bigwig_files_list = utils.load_bigwigs(input_dirs)
    num_bigwigs = len(bigwig_names)
    '''
    bg_list_file = open('%s/bigwig_list.txt'%(input_dirs[0]),'r').readlines()
    big_wig_list = [item.strip() for item in bg_list_file]
    '''
    chip_name_file = np.loadtxt(input_dir + '/chip.txt',dtype=str)
    big_wig_list = [item.split('_')[0]+'.bw' for item in chip_name_file[:,0]]
    model_predicts_list = []
    test_predict_all_TFs = []
    bb = []
    test_predict_all_TFs2 = []
    test_predict_all_TFs3 = []
    bb2 = []
    bb3 = []
    for num_pos,positive_windows in enumerate(positive_windows_list):
        _, meta_names, datagen_bed = utils.load_bed_data_sc(genome, positive_windows, False, False, input_dir, False, big_wig_list,num_pos,chrom=None)
        
        #pdb.set_trace()
        
        print("%d sample...in %d "%(num_pos,len(positive_windows_list)))
        model_predicts,model_tfs = test(output_dir,datagen_bed)
        try:
            os.stat('%s/deepATAC_fea_scbw_moreTFs_weight'%(input_dirs[0]))
        except:
            os.mkdir('%s/deepATAC_fea_scbw_moreTFs_weight'%(input_dirs[0]))
        np.save('%s/deepATAC_fea_scbw_moreTFs_weight/%s_data'%(input_dirs[0],str(num_pos)),model_predicts)
        #model_predicts_list.append(model_predicts) 
        for index in range(model_predicts.shape[0]):
            #test_predict_all_TFs = test_predict_all_TFs + [model_tfs[item] for item in np.argsort(model_predicts[index])[::-1][:1]]
            #test_predict_all_TFs2 = test_predict_all_TFs2 + [model_tfs[item] for item in np.argsort(model_predicts[index])[::-1][:5]]
            test_predict_all_TFs3 = test_predict_all_TFs3 + [model_tfs[item] for item in np.argsort(model_predicts[index])[::-1][:2]]
        #a.append(Counter(test_predict_all_TFs).most_common(12))
        #bb.append([Counter(test_predict_all_TFs)[item]/float(sum(Counter(test_predict_all_TFs).values())) for item in model_tfs])
        #bb2.append([Counter(test_predict_all_TFs2)[item]/float(sum(Counter(test_predict_all_TFs2).values())) for item in model_tfs])
        bb3.append([Counter(test_predict_all_TFs3)[item]/float(sum(Counter(test_predict_all_TFs3).values())) for item in model_tfs])
    #pdb.set_trace()
        #test_predict_all_TFs = []
        #test_predict_all_TFs2 = []
        test_predict_all_TFs3 = []
        
    #pdb.set_trace()
    #bb = np.array(bb)
    #bb2 = np.array(bb2)
    bb3 = np.array(bb3)
    #np.save('%s/count_motif_dis_thres_1_%s'%(input_dirs[0],motif_name),bb)
    #np.save('%s/count_motif_dis_thres_5_%s'%(input_dirs[0],motif_name),bb2)
    np.save('%s/count_motif_scbw_new_dis_thres_2_%s_moreTFs_weight'%(input_dirs[0],motif_name),bb3)
  
if __name__ == '__main__':
    """
    See module-level docstring for a description of the script.
    """
    main()
