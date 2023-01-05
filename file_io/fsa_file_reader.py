
# code taken from https://github.com/bow/abifpy/blob/master/abifpy.py

import datetime
import struct

from sys import version_info

RELEASE = False
__version_info__ = ('1', '0', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


#__all__ = ['Trace']

FILTER_SETS = {
    'G5': {
        '6-FAM':    { 'filter': 'B', 'rgb': (0,0,1) },
        'VIC':      { 'filter': 'G', 'rgb': (0,1,0) },
        'NED':      { 'filter': 'Y', 'rgb': (0.95,0.95,0) },
        'PET':      { 'filter': 'R', 'rgb': (1,0.5,0) },
        'LIZ':      { 'filter': 'O', 'rgb': (1,0,0) }
    }
}

WAVELENGTH = {
    ''
    '6-FAM': 522,
    'VIC': 554,
    'NED': 575,
    'PET': 595,
    'LIZ': 655,
}

# dictionary for deciding which values to extract and contain in self.data
EXTRACT = {
            'TUBE1': 'well',
            'DySN1': 'dye',
            'GTyp1': 'polymer',
            'MODL1': 'model',
            'RUND1': 'run start date',
            'RUND2': 'run finish date',
            'RUND3': 'data collection start date',
            'RUND4': 'data collection finish date',
            'RUNT1': 'run start time',
            'RUNT2': 'run finish time',
            'RUNT3': 'data collection start time',
            'RUNT4': 'data collection finish time',
            'DATA1': 'raw1',
            'DATA2': 'raw2',
            'DATA3': 'raw3',
            'DATA4': 'raw4',
            'PLOC2': 'tracepeaks',
            'FWO_1': 'baseorder',
            'SMap1': 'size map 1',
            'SMap2': 'size map 2',
          }

# dictionary for unpacking tag values
_BYTEFMT = {
            1: 'b',     # byte
            2: 's',     # char
            3: 'H',     # word
            4: 'h',     # short
            5: 'i',     # long
            6: '2i',    # rational, legacy unsupported
            7: 'f',     # float
            8: 'd',     # double
            10: 'h2B',  # date
            11: '4B',   # time
            12: '2i2b', # thumb
            13: 'B',    # bool
            14: '2h',   # point, legacy unsupported
            15: '4h',   # rect, legacy unsupported
            16: '2i',   # vPoint, legacy unsupported
            17: '4i',   # vRect, legacy unsupported
            18: 's',    # pString
            19: 's',    # cString
            20: '2i',   # Tag, legacy unsupported
           }

# header structure
_HEADFMT = '>4sH4sI2H3I'

# directory data structure
_DIRFMT = '>4sI2H4I'

# to handle py3 IO
def py3_get_string(byte):
    if version_info[0] < 3:
        return byte
    else:
        return byte.decode()

def py3_get_byte(string):
    if version_info[0] < 3:
        return string
    else:
        return string.encode()

class TraceDir(object):
    """Class representing directory content."""

    def __init__(self, tag_entry, handle):
        self.tag_name = py3_get_string(tag_entry[0])
        self.tag_num = tag_entry[1]
        self.elem_code = tag_entry[2]
        self.elem_size = tag_entry[3]
        self.elem_num = tag_entry[4]
        self.data_size = tag_entry[5]
        self.data_offset = tag_entry[6]
        self.data_handle = tag_entry[7]
        self.tag_offset = tag_entry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.data_size <= 4:
            self.data_offset = self.tag_offset + 20

        self.tag_data = self._unpack(handle)

    def __repr__(self):
        """Represents data associated with a tag."""
        summary = ['tag_name: {0}'.format(repr(self.tag_name))]
        summary.append('tag_number: {0}'.format(repr(self.tag_num)))
        summary.append('elem_code: {0}'.format(repr(self.elem_code)))
        summary.append('elem_size: {0}'.format(repr(self.elem_size)))
        summary.append('elem_num: {0}'.format(repr(self.elem_num)))
        summary.append('data_size: {0}'.format(repr(self.data_size)))
        summary.append('data_offset: {0}'.format(repr(self.data_offset)))
        summary.append('data_handle: {0}'.format(repr(self.data_handle)))
        summary.append('tag_offset: {0}'.format(repr(self.tag_offset)))
        summary.append('tag_data: {0}'.format(repr(self.tag_data)))

        return '\n'.join(summary)

    def _unpack(self, handle):
        """Returns tag data"""
        if self.elem_code in _BYTEFMT:

            # because ">1s" unpacks differently from ">s"
            num = '' if self.elem_num == 1 else str(self.elem_num)
            fmt = "{0}{1}{2}".format('>', num, _BYTEFMT[self.elem_code])
            start = self.data_offset

            handle.seek(start)
            data = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))

            # no need to use tuple if len(data) == 1
            if self.elem_code not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elem_code == 2:
                return py3_get_string(data)
            elif self.elem_code == 10:
                return datetime.date(*data)
            elif self.elem_code == 11:
                return datetime.time(*data)
            elif self.elem_code == 13:
                return bool(data)
            elif self.elem_code == 18:
                return py3_get_string(data[1:])
            elif self.elem_code == 19:
                return py3_get_string(data[:-1])
            else:
                return data
        else:
            return None

def parse_fsa_header(header, handle):
    """Generator for directory contents."""
    # header structure:
    # file signature, file version, tag name, tag number,
    # element type code, element size, number of elements
    # data size, data offset, handle
    head_elem_size = header[5]
    head_elem_num = header[6]
    head_offset = header[8]
    index = 0

    while index < head_elem_num:
        start = head_offset + index * head_elem_size
        # added directory offset to tuple
        # to handle directories with data size <= 4 bytes
        handle.seek(start)
        dir_entry = struct.unpack(_DIRFMT,
                                  handle.read(struct.calcsize(_DIRFMT))) + (start,)
        index += 1
        yield TraceDir(dir_entry, handle)

def readFSAFile(in_file):
    handle = open(in_file, 'rb')
    #contents = handle.read()
    #print(contents)
    try:
        handle.seek(0)
        results=handle.read(4)
        print(results)
        if not results == py3_get_byte('ABIF'):
            raise IOError('Input is not a valid trace file')
    except IOError:
        handle = None
        raise
    else:
        # header data structure:
        # file type, file, version, tag name, tag number, element type code,
        # element size, number of elements, data size, data offset, handle,
        # file type, file version
        # dictionary for containing file metadata
        print("reading " + in_file)
        data = {}
        # dictionary for containing extracted directory data
        tags = {}
        # values contained in file header
        handle.seek(0)
        header = struct.unpack(_HEADFMT,
                               handle.read(struct.calcsize(_HEADFMT)))
        # file format version
        version = header[1]
        print("FSA file version: " + str(version))

        all_keys = {}
        header_entries = parse_fsa_header(header, handle)

        # print("All Data keys in header: " +  str(keys))

        # build dictionary of data tags and metadata
        for entry in header_entries:
            key = entry.tag_name + str(entry.tag_num)
            tags[key] = entry
            # only extract data from tags we care about
            if key in EXTRACT:
                # e.g. self.data['well'] = 'B6'
                # print("extracting data for key: " + str(key))
                # data[EXTRACT[key]] = get_data(key)
                """Returns data stored in a tag."""
                data[EXTRACT[key]] = tags[key].tag_data

            all_keys[key] = tags[key].tag_data




    handle.close()

    print("Dyes in this file: " )
    dye_to_channel_mapping={}
    channel_numbers=[1,2,3,4]
    for channel_number in channel_numbers:
            channel_name = 'DATA' + str(channel_number)
            wavelength = all_keys['DyeW' + str(channel_number)]
            dye_name = all_keys['DyeN' + str(channel_number)]

            if dye_name =="6-FAM":
                dye_name = "FAM"

            dye_to_channel_mapping[dye_name ]= channel_number
            print(channel_name + ":\t" + dye_name + ", wavelength " +  str(wavelength))

    return dye_to_channel_mapping, all_keys