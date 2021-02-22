#! /usr/bin/python
import sys
import h5py
import os
import fnmatch

reads_group = 'Raw/Reads'
analysis_group = 'Analyses'
# MAX_BASE_NUM = max_num_bases()


def get_fast5s(fast5_dir, is_recursive=True):
    fast5_dir = os.path.abspath(fast5_dir)
    fast5s = []
    if is_recursive:
        for root, dirnames, filenames in os.walk(fast5_dir):
            for filename in fnmatch.filter(filenames, '*.fast5'):
                fast5_path = os.path.join(root, filename)
                fast5s.append(fast5_path)
    else:
        for fast5_name in os.listdir(fast5_dir):
            if fast5_name.endswith('.fast5'):
                fast5_path = '/'.join([fast5_dir, fast5_name])
                fast5s.append(fast5_path)
    return fast5s


def del_group_of_Analyses_in_fast5(filepath, groupnamelist):
    h5file = h5py.File(filepath, mode='r+')
    # grouptmp = h5file[groupname]
    for groupname in groupnamelist:
        grouppath = analysis_group + "/" + groupname
        if grouppath in h5file.keys():
            del h5file[grouppath]
    h5file.close()
    # repack
    repack_bin = "/homeb/nipeng/tools/hdf5-1.12.0/tools/src/h5repack/h5repack"
    parname = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    basename1 = basename + ".repack.hdf5"
    cmd = "{} {} {}".format(repack_bin, parname + "/" + basename, parname + "/" + basename1)
    os.system(cmd)
    cmd = "mv {} {}".format(parname + "/" + basename1, parname + "/" + basename)
    os.system(cmd)


def main(fast5_dir, subg_names):
    subg_name_list = subg_names.split(",")
    fast5s = get_fast5s(fast5_dir, True)
    print("get {} fast5s".format(len(fast5s)))
    cnt = 0
    for fast5file in fast5s:
        del_group_of_Analyses_in_fast5(fast5file, subg_name_list)
        cnt += 1
        print("file {} done..".format(cnt))
    print("done")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
