#! /usr/bin/python
import sys
import h5py
import os
import fnmatch

reads_group = 'Raw/Reads'
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


def get_readid_from_fast5(h5file):
    first_read = list(h5file[reads_group].keys())[0]
    if sys.version_info[0] >= 3:
        read_id = str(h5file['/'.join([reads_group, first_read])].attrs['read_id'], 'utf-8')
    else:
        read_id = str(h5file['/'.join([reads_group, first_read])].attrs['read_id'])
    # print(read_id)
    return read_id


def main(fast5_dir, wfile):
    wf = open(wfile, 'w')
    fast5s = get_fast5s(fast5_dir, True)
    print("get {} fast5s".format(len(fast5s)))
    for fast5file in fast5s:
        h5file = h5py.File(fast5file, mode='r')
        readid = get_readid_from_fast5(h5file)
        wf.write(readid + "\n")
    wf.close()
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
