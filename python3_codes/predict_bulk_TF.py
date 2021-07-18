#!/usr/bin/env python
"""
Script for training model.
python predict_bulk_TF.py-i data/GM12878 -e 5 -oc mul
Use `predict_bulk_TF.py -h` to see an auto-generated description of advanced options.
"""
#import utils_ATAC as utils
#import utils_ATAC_k26 as utils
#import utils_ATAC_add_freq as utils
import utils## use for scFAN training
import numpy as np

# Standard library imports
import sys
import os
import errno
import argparse
import pickle
import pdb

def test(model_dir,datagen_test):
    model_tfs, model_bigwig_names, features, model = utils.load_model(model_dir)
    model_predicts = model.predict_generator(datagen_test, val_samples=1+len(datagen_test)/100, pickle_safe=True,verbose=1)
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
        description="predict bulk TF model.",
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
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-o', '--outputdir', type=str,
                       help='The output directory. Causes error if the directory already exists.')
    group.add_argument('-oc', '--outputdirc', type=str,
                       help='The output directory. The pretrained model directory. Will overwrite if directory already exists.')
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
                print(sys.stderr, ('output directory (%s) already exists '
                                      'but you specified not to clobber it') % output_dir)
                sys.exit(1)
            else:
                print(sys.stderr, ('output directory (%s) already exists '
                                      'so it will be clobbered') % output_dir)

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
        tfs, positive_windows, y_positive, nonnegative_regions_bed = \
            utils.load_chip_multiTask(input_dir)
        num_tfs = len(tfs)
    print('Loading bigWig data')
    bigwig_names, bigwig_files_list = utils.load_bigwigs(input_dirs)
    num_bigwigs = len(bigwig_names)
    if valid_input_dirs:
        valid_bigwig_names, valid_bigwig_files_list = utils.load_bigwigs(valid_input_dirs)
        assert valid_bigwig_names == bigwig_names
    if not singleTask:
        bigwig_files = bigwig_files_list[0]
    if meta:## did not use in scFAN
        print('Loading metadata features')
        
    else:# meta option was not selected, pass empty metadata features to the functions
        meta_list = [[] for bigwig_files in bigwig_files_list]
        if valid_input_dirs:
            valid_meta_list = [[] for bigwig_files in valid_bigwig_files_list]
    
    print('Making features')
    if singleTask:## did not use in scFAN
        print('single Task feature')
    else:
        datagen_train, datagen_valid,datagen_test,data_valid,data_test = \
            utils.make_features_multiTask(positive_windows, y_positive,
            nonnegative_regions_bed, bigwig_files, bigwig_names,
            genome, epochs, valid_chroms, test_chroms)
    #pdb.set_trace()
    print('Building model')
    if len(data_test) == 0:
        model_predicts,model_tfs = test(output_dir,datagen_valid)
        data_test = data_valid
    else:
        model_predicts,model_tfs = test(output_dir,datagen_test)

    ### test and evaluate  users and uncomment this part for evaluation
    '''
    test_label_all = []
    for _,item in enumerate(data_test):
        test_label_all.append(item[-1])
    test_label_all = np.array(test_label_all)
    atac = []
    atacpr = []
    atacrc = []
    #pdb.set_trace()
    for index in range(len(model_tfs)):
        truth = test_label_all[:,index]
        pred = model_predicts[:,index]
        atac.append(roc_auc_score(truth, pred))
        atacpr.append(average_precision_score(truth, pred))
        #atacrc.append(recall_score(truth,pred,average='macro'))
        try:
            precision, recall, thresholds = precision_recall_curve(truth, pred)
            atacrc.append(recall[np.where(precision==0.5)[0][0]])
        except:
            #pdb.set_trace()
            atacrc.append(0)#(recall[np.where((precision>=0.49)&(precision<0.51))[0][0]])
            #continue
    #np.save('draw_plot/%s/more_TFs/auc_value_%s_ATAC_moreTFs_FactorNet'%(motif,motif),atac) #np.save('draw_plot/%s/more_TFs/auc_value_%s_ATAC_moreTFs_FactorNet'%(motif,motif),atac)
    #np.save('draw_plot/%s/more_TFs/aupr_value_%s_ATAC_moreTFs_FactorNet'%(motif,motif),atacpr)
    #np.save('draw_plot/%s/more_TFs/recall_value_%s_ATAC_moreTFs_FactorNet'%(motif,motif),atacrc)
    #pdb.set_trace()
    print("Average AUC ROC (ATAC)", np.mean(atac))
    print("Average AUPR ROC (ATAC)", np.mean(atacpr))
    print("Average Recall ROC (ATAC)", np.mean(atacrc))
    '''
    return model_predicts


if __name__ == '__main__':
    """
    See module-level docstring for a description of the script.
    """
    main()
