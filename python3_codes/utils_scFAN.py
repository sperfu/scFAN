import numpy as np
from pybedtools import BedTool, Interval
import pyfasta
import parmap
import copy
# Standard library imports
import os
import itertools
from data_iter_scFAN import DataIterator
import pybedtools
import pdb
import subprocess

t_flag = subprocess.getstatusoutput('mkdir tmp_scFAN')
pybedtools.set_tempdir('tmp_scFAN')
batch_size = 100

genome_sizes_file = 'resources/hg19.autoX.chrom.sizes'
genome_fasta_file = 'resources/hg19.fa'
blacklist_file = 'resources/blacklist.bed.gz'

genome_window_size = 200
genome_window_step = 50
shift_size = 20

def set_seed(seed):
    np.random.seed(seed)


def chroms_filter(feature, chroms):
    if feature.chrom in chroms:
        return True
    return False


def subset_chroms(chroms, bed):
    result = bed.filter(chroms_filter, chroms).saveas()
    return BedTool(result.fn)


def get_genome_bed():
    genome_sizes_info = np.loadtxt(genome_sizes_file, dtype=str)
    chroms = list(genome_sizes_info[:,0])
    chroms_sizes = list(genome_sizes_info[:,1].astype(int))
    genome_bed = []
    for chrom, chrom_size in zip(chroms, chroms_sizes):
        genome_bed.append(Interval(chrom, 0, chrom_size))
    genome_bed = BedTool(genome_bed)
    return chroms, chroms_sizes, genome_bed


def get_bigwig_rc_order(bigwig_names):
    assert len(set(bigwig_names)) == len(bigwig_names)
    rc_indices = np.arange(len(bigwig_names))
    for ind, bigwig_name in enumerate(bigwig_names):
        if bigwig_name[-4:] == '_fwd':
            bigwig_rc_name = bigwig_name[:-4] + '_rev'
            bigwig_rc_index = bigwig_names.index(bigwig_rc_name)
            rc_indices[bigwig_rc_index] = ind
        if bigwig_name[-4:] == '_rev':
            bigwig_rc_name = bigwig_name[:-4] + '_fwd'
            bigwig_rc_index = bigwig_names.index(bigwig_rc_name)
            rc_indices[bigwig_rc_index] = ind
    return rc_indices


def make_features_multiTask(positive_windows, y_positive, nonnegative_regions_bed, 
                            bigwig_files, bigwig_names, genome, epochs, valid_chroms, test_chroms):
    chroms, chroms_sizes, genome_bed = get_genome_bed()
    train_chroms = chroms
    for chrom in valid_chroms + test_chroms:
        train_chroms.remove(chrom)
    genome_bed_train, genome_bed_valid, genome_bed_test = \
        [subset_chroms(chroms_set, genome_bed) for chroms_set in
         (train_chroms, valid_chroms, test_chroms)]

    positive_windows_train = []
    positive_windows_valid = []
    positive_windows_test = []
    positive_data_train = []
    positive_data_valid = []
    positive_data_test = []
    
    import pdb
    print('Splitting positive windows into training, validation, and testing sets')
    for positive_window, target_array in itertools.izip(positive_windows, y_positive):
        if len(positive_window.chrom) > 8:
            pdb.set_trace()
        chrom = positive_window.chrom
        start = int(positive_window.start)
        stop = int(positive_window.stop)
        if chrom in test_chroms:
            positive_windows_test.append(positive_window)
            positive_data_test.append((chrom, start, stop, shift_size, bigwig_files, [], target_array))
        elif chrom in valid_chroms:
            positive_windows_valid.append(positive_window)
            positive_data_valid.append((chrom, start, stop, shift_size, bigwig_files, [], target_array))
        else:
            positive_windows_train.append(positive_window)
            positive_data_train.append((chrom, start, stop, shift_size, bigwig_files, [], target_array))
    
    positive_windows_train = BedTool(positive_windows_train)
    positive_windows_valid = BedTool(positive_windows_valid)
    positive_windows_test = BedTool(positive_windows_test)

    import pdb
    print('Getting negative training examples')
    negative_windows_train = BedTool.cat(*(epochs*[positive_windows]), postmerge=False)
    #negative_windows_train = BedTool.cat(*(10*[positive_windows]), postmerge=False)
    #pdb.set_trace()
    negative_windows_train = negative_windows_train.shuffle(g=genome_sizes_file,
                                                            incl=genome_bed_train.fn,
                                                            excl=nonnegative_regions_bed.fn,
                                                            noOverlapping=False,
                                                            seed=np.random.randint(-214783648, 2147483647))
                                                            #seed=np.random.randint(-21478364, 21474836))
    print('Getting negative validation examples')
    negative_windows_valid = positive_windows_valid.shuffle(g=genome_sizes_file,
                                                            incl=genome_bed_valid.fn,
                                                            excl=nonnegative_regions_bed.fn,
                                                            noOverlapping=False,
                                                            seed=np.random.randint(-214783648, 2147483647))
                                                            #seed=np.random.randint(-21478364, 21474836))
    print('Getting negative testing examples')
    negative_windows_test = positive_windows_test.shuffle(g=genome_sizes_file,
                                                            incl=genome_bed_test.fn,
                                                            excl=nonnegative_regions_bed.fn,
                                                            noOverlapping=False,
                                                            seed=np.random.randint(-214783648, 2147483647))
                                                            #seed=np.random.randint(-21478364, 21474836))

    # Train
    print('Extracting data from negative training BEDs')
    negative_targets = np.zeros(y_positive.shape[1])
    negative_data_train = [(window.chrom, window.start, window.stop, shift_size, bigwig_files, [], negative_targets)
                           for window in negative_windows_train]

    # Validation
    print('Extracting data from negative validation BEDs')
    negative_data_valid = [(window.chrom, window.start, window.stop, shift_size, bigwig_files, [], negative_targets)
                           for window in negative_windows_valid]
    
    # Test
    print('Extracting data from negative testing BEDs')
    negative_data_test = [(window.chrom, window.start, window.stop, shift_size, bigwig_files, [], negative_targets)
                           for window in negative_windows_test]

    num_positive_train_windows = len(positive_data_train)
    
    data_valid = negative_data_valid + positive_data_valid
    data_test = negative_data_test + positive_data_test

    print('Shuffling training data')
    data_train = []
    for i in xrange(epochs):
        epoch_data = []
        epoch_data.extend(positive_data_train)
        epoch_data.extend(negative_data_train[i*num_positive_train_windows:(i+1)*num_positive_train_windows])
        np.random.shuffle(epoch_data)
        data_train.extend(epoch_data)

    print('Generating data iterators')
    bigwig_rc_order = get_bigwig_rc_order(bigwig_names)
    datagen_train = DataIterator(data_train, genome, batch_size, L, bigwig_rc_order)
    datagen_valid = DataIterator(data_valid, genome, batch_size, L, bigwig_rc_order)
    datagen_test = DataIterator(data_test, genome, batch_size, L, bigwig_rc_order)

    print(len(datagen_train), 'training samples')
    print(len(datagen_valid), 'validation samples')
    print(len(datagen_test), 'test samples')
    return datagen_train, datagen_valid, datagen_test, data_test


def data_to_bed(data):
    intervals = []
    for datum in data:
        chrom = datum[0]
        start = datum[1]
        stop = datum[2]
        intervals.append(Interval(chrom, start, stop))
    return BedTool(intervals)


def extract_data_from_bed(args, shift, label, gencode):
    peaks = args[0]
    bigwig_files = args[1]
    meta = args[2]

    data = []    

    if gencode:
        cpg_bed = BedTool('resources/cpgisland.bed.gz')
        cds_bed = BedTool('resources/wgEncodeGencodeBasicV19.cds.merged.bed.gz')
        intron_bed = BedTool('resources/wgEncodeGencodeBasicV19.intron.merged.bed.gz')
        promoter_bed = BedTool('resources/wgEncodeGencodeBasicV19.promoter.merged.bed.gz')
        utr5_bed = BedTool('resources/wgEncodeGencodeBasicV19.utr5.merged.bed.gz')
        utr3_bed = BedTool('resources/wgEncodeGencodeBasicV19.utr3.merged.bed.gz')

        peaks_cpg_bedgraph = peaks.intersect(cpg_bed, wa=True, c=True)
        peaks_cds_bedgraph = peaks.intersect(cds_bed, wa=True, c=True)
        peaks_intron_bedgraph = peaks.intersect(intron_bed, wa=True, c=True)
        peaks_promoter_bedgraph = peaks.intersect(promoter_bed, wa=True, c=True)
        peaks_utr5_bedgraph = peaks.intersect(utr5_bed, wa=True, c=True)
        peaks_utr3_bedgraph = peaks.intersect(utr3_bed, wa=True, c=True)

        for cpg, cds, intron, promoter, utr5, utr3 in itertools.izip(peaks_cpg_bedgraph,peaks_cds_bedgraph,peaks_intron_bedgraph,peaks_promoter_bedgraph,peaks_utr5_bedgraph,peaks_utr3_bedgraph):
            chrom = cpg.chrom
            peak_start = cpg.start
            peak_stop = cpg.stop
            peak_mid = (peak_start + peak_stop)/2
            start = peak_mid - genome_window_size/2
            stop = peak_mid + genome_window_size/2
            if shift:
                shift_size = peak_stop - start - 75 - 1
            else:
                shift_size = 0
            gencode = np.array([cpg.count, cds.count, intron.count, promoter.count, utr5.count, utr3.count], dtype=bool)
            meta_gencode = np.append(meta, gencode)
            data.append((chrom, start, stop, shift_size, bigwig_files, meta_gencode, label))
    else:
        for peak in peaks:
            chrom = peak.chrom
            peak_start = peak.start
            peak_stop = peak.stop
            peak_mid = (peak_start + peak_stop)/2
            start = peak_mid - genome_window_size/2
            stop = peak_mid + genome_window_size/2
            if shift:
                shift_size = peak_stop - start - 75 - 1
            else:
                shift_size = 0
            data.append((chrom, start, stop, shift_size, bigwig_files, meta, label))

    return data


def valid_test_split_wrapper(bed, valid_chroms, test_chroms):
    bed_train = []
    bed_valid = []
    bed_test = []
    for interval in bed:
        chrom = interval.chrom
        start = interval.start
        stop = interval.stop
        if chrom in test_chroms:
            bed_test.append(interval)
        elif chrom in valid_chroms:
            bed_valid.append(interval)
        else:
            bed_train.append(interval)
    bed_train = BedTool(bed_train)
    bed_valid = BedTool(bed_valid)
    bed_test = BedTool(bed_test)
    return bed_train, bed_valid, bed_test


def negative_shuffle_wrapper(args, include_bed, num_copies, noOverlapping):
    positive_windows = args[0]
    nonnegative_regions_bed = args[1]
    bigwig_files = args[2]
    randomseed = args[3]
    if num_copies > 1:
        positive_windows = BedTool.cat(*(num_copies * [positive_windows]), postmerge=False)
    negative_windows = positive_windows.shuffle(g=genome_sizes_file,
                                                incl=include_bed.fn,
                                                excl=nonnegative_regions_bed.fn,
                                                noOverlapping=noOverlapping,
                                                seed=randomseed)
    return negative_windows

def get_onehot_chrom(chrom): 
    fasta = pyfasta.Fasta(genome_fasta_file)
    chr_str = str(fasta[chrom]).upper()
    d = np.array(['A','C','G','T'], dtype='|S1')
    y = np.fromstring(chr_str, dtype='|S1')[:, np.newaxis] == d
    return y


def load_genome():
    chroms = list(np.loadtxt(genome_sizes_file, usecols=[0], dtype=str))
    onehot_chroms = parmap.map(get_onehot_chrom, chroms)
    genome_dict = dict(zip(chroms, onehot_chroms))
    return genome_dict


def intersect_count(chip_bed, windows_file):
    windows = BedTool(windows_file)
    chip_bedgraph = windows.intersect(chip_bed, wa=True, c=True, f=1.0*(genome_window_size/2+1)/genome_window_size, sorted=True)
    bed_counts = [i.count for i in chip_bedgraph]
    return bed_counts


def make_blacklist():
    blacklist = BedTool(blacklist_file)
    blacklist = blacklist.slop(g=genome_sizes_file, b=L)
    # Add ends of the chromosomes to the blacklist
    genome_sizes_info = np.loadtxt(genome_sizes_file, dtype=str)
    chroms = list(genome_sizes_info[:,0])
    chroms_sizes = list(genome_sizes_info[:,1].astype(int))
    blacklist2 = []
    for chrom, size in zip(chroms, chroms_sizes):
        blacklist2.append(Interval(chrom, 0, L))
        blacklist2.append(Interval(chrom, size - L, size))
    blacklist2 = BedTool(blacklist2)
    blacklist = blacklist.cat(blacklist2)
    return blacklist

def get_chip_beds_multiple(input_dir,process_batch):
    chip_info_file = input_dir + '/chip.txt'
    chip_info = np.loadtxt(chip_info_file, dtype=str)
    if len(chip_info.shape) == 1:
        chip_info = np.reshape(chip_info, (-1,len(chip_info)))
    tfs = list(chip_info[:, 1])
    chip_bed_files = [input_dir + '/' + i for i in chip_info[:,0]]
    chip_beds = [BedTool(chip_bed_file) for chip_bed_file in chip_bed_files]
    if process_batch:
        batch_name_list = list(np.unique(chip_info[:,1]))
        batch_list_all_dict = {}
        exchange_dict = {}
        for index,item in enumerate(batch_name_list):
            batch_tmp = [chip_beds[i] for i in list(np.where(chip_info[:,1]==item)[0])]
            batch_0 = batch_tmp[0]
            batch_tmp = batch_tmp[1:]
            print('concatenate batch bedfiles for batch %d...'%(index))
            batch_list = batch_0.cat(*batch_tmp,postmerge=False)
            if item not in batch_list_all_dict.keys():
                batch_list_all_dict[item] = batch_list
                batch_name_list_tmp = copy.deepcopy(batch_name_list)
                batch_name_list_tmp.remove(item)
                if len(batch_name_list_tmp) > 1:
                    exchange_dict[item] = batch_name_list_tmp
                else:
                    exchange_dict[item] = batch_name_list_tmp[0]
            else:
                print("Error!!!")
    else:
        print('No need process batch,continue...')
    print('Sorting BED files')
    chip_beds = [chip_bed.sort() for chip_bed in chip_beds]
    merged_chip_bed_list = []
    for item in chip_beds:
        if 1 > 1:
            merged_chip_bed = BedTool.cat(*item)
        else:
            merged_chip_bed = item
        merged_chip_bed_list.append(merged_chip_bed)
    if process_batch:
        return tfs, chip_beds, merged_chip_bed_list,batch_list_all_dict,chip_info[:,1],exchange_dict
    else:
        return tfs, chip_beds, merged_chip_bed_list


def get_chip_beds(input_dir):
    chip_info_file = input_dir + '/chip.txt'
    chip_info = np.loadtxt(chip_info_file, dtype=str)
    if len(chip_info.shape) == 1:
        chip_info = np.reshape(chip_info, (-1,len(chip_info)))
    tfs = list(chip_info[:, 1])
    chip_bed_files = [input_dir + '/' + i for i in chip_info[:,0]]
    pdb.set_trace()
    chip_beds = [BedTool(chip_bed_file) for chip_bed_file in chip_bed_files]
    print('Sorting BED files')
    chip_beds = [chip_bed.sort() for chip_bed in chip_beds]
    if len(chip_beds) > 1:
        merged_chip_bed = BedTool.cat(*chip_beds)
    else:
        merged_chip_bed = chip_beds[0]
    return tfs, chip_beds, merged_chip_bed

def generate_dict_for_batch_cat(batch_list_all_dict,batch_info,exchange_dict):
    batch_tmp_all_dict = {}
    for item in np.unique(batch_info):
        print('Process batch %s...'%(item))
        #pdb.set_trace()
        if item not in batch_tmp_all_dict.keys():
            batch_list_tmp = exchange_dict[item]
            batch_tmp_0 = batch_list_all_dict[batch_list_tmp[0]]
            for tmp in batch_list_tmp[1:]:
                batch_tmp_all = batch_tmp_0.cat(batch_list_all_dict[tmp],postmerge=False)
                batch_tmp_0 = batch_tmp_all
            batch_tmp_all_dict[item] = batch_tmp_0
        else:
            print("duplicate...")
    print('Done merge batch...')
    return batch_tmp_all_dict

def load_chip_multiTask_multiple(input_dir,process_batch):
    if process_batch:
        tfs, chip_beds_list, merged_chip_bed_list, batch_list_all_dict, batch_info, exchange_dict = get_chip_beds_multiple(input_dir,process_batch)
    else:
        tfs, chip_beds_list, merged_chip_bed_list = get_chip_beds_multiple(input_dir,process_batch)
    print('Removing peaks outside of X chromosome and autosomes')
    chroms, chroms_sizes, genome_bed = get_genome_bed()
    blacklist = make_blacklist()
    print('Prepare for all the merge batch files')
    if process_batch and len(exchange_dict) > 2:
        batch_tmp_all_dict = generate_dict_for_batch_cat(batch_list_all_dict,batch_info,exchange_dict)
    positive_windows_list = []
    print('Windowing genome')
    genome_windows = BedTool().window_maker(g=genome_sizes_file, w=genome_window_size,
                                            s=genome_window_step)
    for index,merged_chip_bed in enumerate(merged_chip_bed_list):
        if process_batch:
            print('processing batch peaks...')
            #pdb.set_trace()
            if len(exchange_dict) == 2:
                #merged_chip_bed = merged_chip_bed.intersect(batch_list_all_dict[exchange_dict[batch_info[index]]],wa=True,u=True)
                merged_chip_bed = merged_chip_bed.intersect(batch_list_all_dict[exchange_dict[batch_info[index]]])
                merged_chip_bed = merged_chip_bed.sort().merge()
                merged_chip_bed = merged_chip_bed.intersect(genome_bed, u=True, sorted=True)
            else:
                merged_chip_bed = merged_chip_bed.intersect(batch_tmp_all_dict[batch_info[index]])
                merged_chip_bed = merged_chip_bed.sort().merge()
                merged_chip_bed = merged_chip_bed.intersect(genome_bed, u=True, sorted=True)
        else:
            merged_chip_bed = merged_chip_bed.intersect(genome_bed, u=True, sorted=True)

        print('Dataset %d/%d...'%(index+1,len(merged_chip_bed_list)))

        print('Extracting windows that overlap at least one ChIP interval')
        positive_windows = genome_windows.intersect(merged_chip_bed, u=True, f=1.0*(genome_window_size/2+1)/genome_window_size, sorted=True)

        # Exclude all windows that overlap a blacklisted region
        
        print('Removing windows that overlap a blacklisted region')
        positive_windows = positive_windows.intersect(blacklist, wa=True, v=True, sorted=True)

        # Binary binding target matrix of all positive windows
        #pdb.set_trace()
        # Generate targets
        # Later we want to gather negative windows from the genome that do not overlap
        # with a blacklisted or ChIP region
        positive_windows_list.append(positive_windows)
    return tfs, positive_windows_list


def nonnegative_wrapper(a, bl_file):
    bl = BedTool(bl_file)
    a_slop = a.slop(g=genome_sizes_file, b=genome_window_size)
    return bl.cat(a_slop).fn


def get_chip_bed(input_dir, tf, bl_file):
    blacklist = BedTool(bl_file)
    chip_info_file = input_dir + '/chip.txt'
    chip_info = np.loadtxt(chip_info_file, dtype=str)
    if len(chip_info.shape) == 1:
        chip_info = np.reshape(chip_info, (-1,len(chip_info)))
    tfs = list(chip_info[:, 1])
    assert tf in tfs
    tf_index = tfs.index(tf)
    chip_bed_file = input_dir + '/' + chip_info[tf_index, 0]
    chip_bed = BedTool(chip_bed_file)
    chip_bed = chip_bed.sort()
    #Remove any peaks not in autosomes or X chromosome
    chroms, chroms_sizes, genome_bed = get_genome_bed()
    chip_bed = chip_bed.intersect(genome_bed, u=True, sorted=True)
    #Remove any peaks in blacklist regions
    chip_bed = chip_bed.intersect(blacklist, wa=True, v=True, sorted=True)
    if chip_info.shape[1] == 3:
        relaxed_bed_file = input_dir + '/' + chip_info[tf_index, 2]
        relaxed_bed = BedTool(relaxed_bed_file)
        relaxed_bed = relaxed_bed.sort()
    else:
        relaxed_bed = chip_bed
    return chip_bed, relaxed_bed

def load_bigwigs_sc(input_dirs,num_pos,bigwig_lists,input_scATAC_dir):
    bigwig_names = None
    bigwig_files_list = []
    for input_dir in input_dirs:
        input_bigwig_info_file = input_dir + '/bigwig.txt'
        if not os.path.isfile(input_bigwig_info_file):
            input_bigwig_names = []
            input_bigwig_files_list = []
        else:
            input_bigwig_info = np.loadtxt(input_bigwig_info_file, dtype=str)
            if len(input_bigwig_info.shape) == 1:
                input_bigwig_info = np.reshape(input_bigwig_info, (-1,2))
            input_bigwig_names = list(input_bigwig_info[:, 1])
            ### changed here for single cell
            input_bigwig_files = [input_dir + '/' + i for i in input_bigwig_info[:,0][:-1]] ## previous_version
            #input_bigwig_files = [input_dir + '/' + i for i in input_bigwig_info[:,0]]

            input_bigwig_files.append(input_scATAC_dir[0]+'/'+bigwig_lists[num_pos])
            #input_bigwig_files.append('/data2/fly/scFAN_data/bigwig_all/'+bigwig_lists[num_pos])
            print(input_bigwig_files)
            #pdb.set_trace()
        if bigwig_names is None:
            bigwig_names = input_bigwig_names
        else:
            assert bigwig_names == input_bigwig_names
        bigwig_files_list.append(input_bigwig_files)
    return bigwig_names, bigwig_files_list

def load_bigwigs(input_dirs):
    bigwig_names = None
    bigwig_files_list = []
    for input_dir in input_dirs:
        input_bigwig_info_file = input_dir + '/bigwig.txt'
        if not os.path.isfile(input_bigwig_info_file):
            input_bigwig_names = []
            input_bigwig_files_list = []
        else:
            input_bigwig_info = np.loadtxt(input_bigwig_info_file, dtype=str)
            if len(input_bigwig_info.shape) == 1:
                input_bigwig_info = np.reshape(input_bigwig_info, (-1,2))
            input_bigwig_names = list(input_bigwig_info[:, 1])
            input_bigwig_files = [input_dir + '/' + i for i in input_bigwig_info[:,0]]
        if bigwig_names is None:
            bigwig_names = input_bigwig_names
        else:
            assert bigwig_names == input_bigwig_names
        bigwig_files_list.append(input_bigwig_files)
    return bigwig_names, bigwig_files_list


def get_output(input_layer, hidden_layers):
    output = input_layer
    for hidden_layer in hidden_layers:
        output = hidden_layer(output)
    return output

def scFANet(out,num_recurrent,num_bws):
    from keras.models import Sequential
    from keras.layers import Flatten, Dense, Dropout, Merge
    from keras.layers import Reshape
    from keras.layers.core import Activation
    from keras.layers.normalization import BatchNormalization
    from keras.layers.convolutional import UpSampling2D
    from keras.layers.convolutional import Conv2D, MaxPooling2D
    from keras.layers.pooling import GlobalMaxPooling2D
    from keras.layers.core import Flatten
    from keras.optimizers import SGD
    from keras import regularizers
    from keras.callbacks import ModelCheckpoint
    from keras.layers import Concatenate
    from keras.models import Model
    from keras.layers import Input
    
    forward_input = Input(shape=(1, L, 4 + num_bws*2,))
    #reverse_input = Input(shape=(L, 4 + num_bws,))
    nkernels = [160,240,480]
    in_size = (1,1000,6)
    l2_lam = 5e-07
    l1_lam = 1e-08
    if num_recurrent > 0:
        hidden_layers = [
            Conv2D(nkernels[0], kernel_size=(1,8), strides=(1,1), padding='same', input_shape=in_size, kernel_regularizer=regularizers.l2(l2_lam)),
            BatchNormalization(),
            Activation('relu'),
            MaxPooling2D(pool_size=(1,4), strides=(1,4)),
            Dropout(0.2),
            
            Conv2D(nkernels[1], kernel_size=(1,8), strides=(1,1), padding='same', kernel_regularizer=regularizers.l2(l2_lam)),
            BatchNormalization(),
            Activation('relu'),
            MaxPooling2D(pool_size=(1,4), strides=(1,4)),
            Dropout(0.2),

            Conv2D(nkernels[1], kernel_size=(1,8), strides=(1,1), padding='same', kernel_regularizer=regularizers.l2(l2_lam)),
            BatchNormalization(),
            Activation('relu'),
            Dropout(0.5),

            Flatten(),
            Dense(919, kernel_regularizer=regularizers.l1(l1_lam)),
            Activation('relu'),
            Dense(out, kernel_regularizer=regularizers.l1(l1_lam)),
            Activation('sigmoid')
        ]
    else:
        pdb.set_trace()
    forward_output = get_output(forward_input, hidden_layers)     
    #reverse_output = get_output(reverse_input, hidden_layers)
    #output = merge([forward_output, reverse_output], mode='ave')
    output = forward_output
    #model = Model(input=[forward_input, reverse_input], output=output)
    model = Model(input=forward_input, output=output)

    return model

def make_model(num_tfs, num_bws, num_motifs, num_recurrent, num_dense, dropout_rate):
    from keras import backend as K
    from keras.models import Model
    from keras.layers import Dense, Dropout, Activation, Flatten, Layer, merge, Input
    from keras.layers.convolutional import Convolution1D, MaxPooling1D
    from keras.layers.pooling import GlobalMaxPooling1D
    from keras.layers.recurrent import LSTM
    from keras.layers.wrappers import Bidirectional, TimeDistributed
    '''
    import tensorflow as tf
    config = tf.ConfigProto(device_count={'gpu':1})
    config.gpu_options.allow_growth=True
    session = tf.Session(config=config)
    '''
    forward_input = Input(shape=(L, 4 + num_bws,))
    reverse_input = Input(shape=(L, 4 + num_bws,))
    if num_recurrent < 0:
        hidden_layers = [
            Convolution1D(input_dim=4 + num_bws, nb_filter=num_motifs,
                          filter_length=w, border_mode='valid', activation='relu',
                          subsample_length=1),
            Dropout(0.1),
            TimeDistributed(Dense(num_motifs, activation='relu')),
            GlobalMaxPooling1D(),
            Dropout(dropout_rate),
            Dense(num_dense, activation='relu'),
            Dropout(dropout_rate),
            Dense(num_tfs, activation='sigmoid')
        ]
    elif num_recurrent == 0:
        hidden_layers = [
            Convolution1D(input_dim=4 + num_bws, nb_filter=num_motifs,
                          filter_length=w, border_mode='valid', activation='relu',
                          subsample_length=1),
            Dropout(0.1),
            TimeDistributed(Dense(num_motifs, activation='relu')),
            MaxPooling1D(pool_length=w2, stride=w2),
            Dropout(dropout_rate),
            Flatten(),
            Dense(num_dense, activation='relu'),
            Dropout(dropout_rate),
            Dense(num_tfs, activation='sigmoid')
        ]
    else:
        hidden_layers = [
            Convolution1D(input_dim=4 + num_bws, nb_filter=num_motifs,
                          filter_length=w, border_mode='valid', activation='relu',
                          subsample_length=1),
            Dropout(0.1),
            TimeDistributed(Dense(num_motifs, activation='relu')),
            MaxPooling1D(pool_length=w2, stride=w2),
            Bidirectional(LSTM(num_recurrent, dropout_W=0.1, dropout_U=0.1, return_sequences=True)),
            Dropout(dropout_rate),
            Flatten(),
            Dense(num_dense, activation='relu'),
            Dropout(dropout_rate),
            Dense(num_tfs, activation='sigmoid')
        ]
    forward_output = get_output(forward_input, hidden_layers)     
    reverse_output = get_output(reverse_input, hidden_layers)
    output = merge([forward_output, reverse_output], mode='ave')
    model = Model(input=[forward_input, reverse_input], output=output)

    return model

def load_model(modeldir):
    tfs_file = modeldir + '/chip.txt'
    tfs = np.loadtxt(tfs_file, dtype=str)
    if len(tfs.shape) == 0:
        tfs = [str(tfs)]
    else:
        tfs = list(tfs)
    bigwig_names_file = modeldir + '/bigwig.txt'
    if not os.path.isfile(bigwig_names_file):
        bigwig_names = []
    else:
        bigwig_names = np.loadtxt(bigwig_names_file, dtype=str)
        if len(bigwig_names.shape) == 0:
            bigwig_names = [str(bigwig_names)]
        else:
            bigwig_names = list(bigwig_names)
    features_file = modeldir + '/feature.txt'
    assert os.path.isfile(features_file)
    features = np.loadtxt(features_file, dtype=str)
    if len(features.shape) == 0:
        features = [str(features)]
    else:
        features = list(features)
    from keras.models import model_from_json
    model_json_file = open(modeldir + '/model.json', 'r')
    model_json = model_json_file.read()
    model = model_from_json(model_json)
    model.load_weights(modeldir + '/best_model.hdf5')
    return tfs, bigwig_names, features, model


def load_bed_data_sc(genome, positive_windows, use_meta, use_gencode, input_dir, is_sorted, big_wig_list, num_pos, input_scATAC_dir, chrom=None):
    bed_filtered = positive_windows
    print('Generating test data iterator')
    bigwig_names, bigwig_files_list  = load_bigwigs_sc([input_dir],num_pos,big_wig_list,input_scATAC_dir)
    bigwig_files = bigwig_files_list[0]
    if use_meta:
        meta_names, meta_list = load_meta([input_dir])
        meta = meta_list[0]
    else:
        meta = []
        meta_names = None

    shift = 0

    if use_gencode:
        cpg_bed = BedTool('resources/cpgisland.bed.gz')
        cds_bed = BedTool('resources/wgEncodeGencodeBasicV19.cds.merged.bed.gz')
        intron_bed = BedTool('resources/wgEncodeGencodeBasicV19.intron.merged.bed.gz')
        promoter_bed = BedTool('resources/wgEncodeGencodeBasicV19.promoter.merged.bed.gz')
        utr5_bed = BedTool('resources/wgEncodeGencodeBasicV19.utr5.merged.bed.gz')
        utr3_bed = BedTool('resources/wgEncodeGencodeBasicV19.utr3.merged.bed.gz')

        peaks_cpg_bedgraph = bed_filtered.intersect(cpg_bed, wa=True, c=True)
        peaks_cds_bedgraph = bed_filtered.intersect(cds_bed, wa=True, c=True)
        peaks_intron_bedgraph = bed_filtered.intersect(intron_bed, wa=True, c=True)
        peaks_promoter_bedgraph = bed_filtered.intersect(promoter_bed, wa=True, c=True)
        peaks_utr5_bedgraph = bed_filtered.intersect(utr5_bed, wa=True, c=True)
        peaks_utr3_bedgraph = bed_filtered.intersect(utr3_bed, wa=True, c=True)

        data_bed = [(window.chrom, window.start, window.stop, 0, bigwig_files, np.append(meta, np.array([cpg.count, cds.count, intron.count, promoter.count, utr5.count, utr3.count], dtype=bool)))
                    for window, cpg, cds, intron, promoter, utr5, utr3 in
                    itertools.izip(bed_filtered, peaks_cpg_bedgraph,peaks_cds_bedgraph,peaks_intron_bedgraph,peaks_promoter_bedgraph,peaks_utr5_bedgraph,peaks_utr3_bedgraph)]
    else:
        data_bed = [(window.chrom, window.start, window.stop, shift, bigwig_files, meta)
                    for window in bed_filtered]
    #from data_iter import DataIterator
    from data_iter_scFAN import DataIterator
    #pdb.set_trace()
    #tmpp = bed_filtered.saveas('/data1/fly/FactorNet/draw_plot/heatmap_plot/merge_heatmap/H1_bed/%d_bed.bed'%(num_pos))
    bigwig_rc_order = get_bigwig_rc_order(bigwig_names)
    datagen_bed = DataIterator(data_bed, genome, 100, L, bigwig_rc_order, shuffle=False)
    return bigwig_names, datagen_bed
