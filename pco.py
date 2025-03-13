import numpy as np
import cv2

"""
Decode PCO image timestamps from binary-coded decimal (see p94 of
"pco_camera_control_commands_105.pdf"). In this version of 'packed BCD' each
pixel contains 2 digits of information in a single byte (8bits). The lower
and upper nibbles (2 x 4bits) encode the numbers 0-9 which are then combined
to give a value in the range 0-99.
"""

def decode_timestamp(image_path):
    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED).astype('uint16')
    # print(image[0,:14])
    # print(image.dtype)    
    # assert len(image.shape) == 2 and image.dtype == 'uint16'
    bcd_px = image[0,:14]                      # get BCD pixels
    lower_nibbles =  bcd_px & 0b00001111        # get lower nibbles
    upper_nibbles = (bcd_px & 0b11110000) >> 4  # get upper nibbles and shift
    dec_px = 10 * upper_nibbles + lower_nibbles # convert to decimal
    # print(dec_px)
    timestamp = {}
    timestamp['image_num'] = np.sum(
        dec_px[:4] * np.array((1e6, 1e4, 1e2, 1)), dtype='uint16')
    timestamp['DD'] = dec_px[7].astype('uint32')
    timestamp['MM'] = dec_px[6].astype('uint32')    
    timestamp['YYYY'] = np.sum(
        dec_px[4:6] * np.array((1e2, 1)), dtype='uint32')
    timestamp['h'] = dec_px[8].astype('uint32')
    timestamp['min'] = dec_px[9].astype('uint32')
    timestamp['s'] = dec_px[10].astype('uint32')
    timestamp['us'] = np.sum(
        dec_px[11:14] * np.array((1e4, 1e2, 1), dtype='uint64'))
    timestamp['time_us'] = np.sum(              # total us on a given day
        dec_px[8:14] * np.array((36e8, 60e6, 1e6, 1e4, 1e2, 1)), dtype='uint64')
    return timestamp



# print(decode_timestamp("E:\\shake_table_data\\N4\\14hz_hopperflow\\60deg\\10cm\\1\\1_00001.tif"))
# print(decode_timestamp("E:\\shake_table_data\\N4\\14hz_hopperflow\\60deg\\10cm\\1\\1_00002.tif"))