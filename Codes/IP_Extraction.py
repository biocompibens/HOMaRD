import tifffile
import numpy as np
import pandas as pd
import skimage
from skimage import transform


def Split_channels_Illu(img_Tiff, cycle_number, path_illu):  # Read and split the TIFF image

    yoyo_channel = []
    rep_channel = []
    nick_channel = []

    original_images_1140 = tifffile.imread(img_Tiff)  # read the TIFF image
    nb_channels = (len(original_images_1140) / 1140)  # Count the number of channels in the TIFF image

    all_illu = np.load(path_illu)
    illu_yoyo = all_illu[0][cycle_number - 1]
    illu_rep = all_illu[1][cycle_number - 1]
    illu_nick = all_illu[2][cycle_number - 1]

    # Flip the original image horizontally and vertically
    if nb_channels == 2:  # If the number of channel is 2
        for i in range(0, len(original_images_1140), 2):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])) - illu_yoyo)
            nick_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])) - illu_nick)
        print "FLIP Splitting Done 2 channels"
        return yoyo_channel, nick_channel
    if nb_channels == 3:  # If the number of channel is 3
        for i in range(0, len(original_images_1140), 3):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])) - illu_yoyo)
            rep_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])) - illu_rep)
            nick_channel.append(np.flipud(np.fliplr(original_images_1140[i + 2])) - illu_nick)
        print "FLIP Splitting Done 3 channels "
        return yoyo_channel, rep_channel, nick_channel


def IP_extraction_subgroup(row, transform_matrix_rep, transform_matrix_nick, yoyo_TIFF, nick_TIFF, rep_TIFF,
                           nb_channel):
    molecule_IP_Yoyo = []
    molecule_IP_Rep = []
    molecule_IP_Nick = []

    list_avg_mol_valueYOYO = []
    list_avg_mol_valueREP = []
    list_avg_mol_valueNICK = []

    transform_matrix_translation_rep = transform.SimilarityTransform(
        translation=[transform_matrix_rep[6], transform_matrix_rep[7]],
        scale=transform_matrix_rep[0])
    transform_matrix_translation_nick = transform.SimilarityTransform(
        translation=[transform_matrix_nick[6], transform_matrix_nick[7]],
        scale=transform_matrix_nick[0])

    if nb_channel == 3:
        if row['RowStart'] == row['RowEnd']:  # For molecules that are in the same FOV

            image_frame = int(row['FovStart'] - 1)  # Frame where the molecule starts
            ori_yoyo_img = yoyo_TIFF[image_frame]
            ori_nick_img = nick_TIFF[image_frame]
            ori_rep_img = rep_TIFF[image_frame]

            # YOYO padding
            image_IP_Yoyo = np.pad(ori_yoyo_img, ((2, 2), (2, 2)), 'constant')
            image_IP_Rep0 = skimage.transform.warp(ori_rep_img, transform_matrix_translation_rep.inverse,
                                                   preserve_range=True, order=3)
            image_IP_Rep = np.pad(image_IP_Rep0, ((2, 2), (2, 2)), 'constant')
            image_IP_Nick0 = skimage.transform.warp(ori_nick_img, transform_matrix_translation_nick.inverse,
                                                    preserve_range=True, order=3)
            image_IP_Nick = np.pad(image_IP_Nick0, ((2, 2), (2, 2)), 'constant')

            begin = int(round(row['YStart'] + 2))  # Y Start position (Round the value and then convert it to int)
            end = int(round(
                row['YEnd'] + 2))  # Y end position (+1 because of the ranges) +2 padding +1 include the last point
            X_position_init = int((round(row['XStart'] + 2)))  # X start

            y_pos_iter = 0
            for Y_pos in range(begin, end + 1):  # +1 because of the padding done previously

                # Blue Red & Green channel extraction IP coord.
                sum_intensity_yoyo = 0
                sum_intensity_rep = 0
                sum_intensity_nick = 0

                X_position = int((X_position_init + (y_pos_iter / row['MoleculeSlope'])))
                y_pos_iter += 1

                for ind_move in [-1, 0, 1]:
                    mol_pixel_yoyo = int(image_IP_Yoyo[Y_pos][X_position + ind_move])
                    sum_intensity_yoyo += mol_pixel_yoyo

                    mol_pixel_rep = int(image_IP_Rep[Y_pos][X_position + ind_move])
                    sum_intensity_rep += mol_pixel_rep

                    mol_pixel_nick = int(image_IP_Nick[Y_pos][X_position + ind_move])
                    sum_intensity_nick += mol_pixel_nick

                avg_mol_pixel_YOYO = sum_intensity_yoyo / 3
                list_avg_mol_valueYOYO.append(avg_mol_pixel_YOYO)

                avg_mol_pixel_REP = sum_intensity_rep / 3
                list_avg_mol_valueREP.append(avg_mol_pixel_REP)

                avg_mol_pixel_NICK = sum_intensity_nick / 3
                list_avg_mol_valueNICK.append(avg_mol_pixel_NICK)

            molecule_IP_Yoyo.append(list_avg_mol_valueYOYO)
            molecule_IP_Nick.append(list_avg_mol_valueNICK)
            molecule_IP_Rep.append(list_avg_mol_valueREP)
        else:  # Molecules overlapping on more than 1 frame
            nbr_overlapping_frame = row['RowEnd'] - row['RowStart']  # Get on how many frames the molecule overlap

            if nbr_overlapping_frame != 1:  # If the number of overlapped frame is greater than 1 (overlap on more than 2 frames)
                molecule_IP_Yoyo.append(0)
                molecule_IP_Rep.append(0)
                molecule_IP_Nick.append(0)

            else:  # If the number of overlapped frame is equal to 1
                molecule_IP_Yoyo.append(0)
                molecule_IP_Rep.append(0)
                molecule_IP_Nick.append(0)
    return pd.Series([row['MoleculeId'], molecule_IP_Yoyo[0], molecule_IP_Rep[0], molecule_IP_Nick[0]])


def IP_extraction(subgroup, transform_matrix_rep, transform_matrix_nick, path_illu, path_no_extension):
    cycle_number = int(subgroup['Scan'].get_values()[0])  # Get the cycle number to open the correct TIFF image
    img_path_init = path_no_extension
    # Split the original TIFF image
    splitted_image = Split_channels_Illu(img_path_init, cycle_number, path_illu)

    if len(splitted_image) == 3:
        yoyo_TIFF = splitted_image[0]
        rep_TIFF = splitted_image[1]
        nick_TIFF = splitted_image[2]
        nb_channel = len(splitted_image)
    else:
        yoyo_TIFF = splitted_image[0]
        nick_TIFF = splitted_image[1]
        rep_TIFF = [0]
        nb_channel = len(splitted_image)

    # Extract each IP for each line of the subgroup
    IP_sub = subgroup.apply(IP_extraction_subgroup, args=(transform_matrix_rep,
                                                          transform_matrix_nick,
                                                          yoyo_TIFF,
                                                          nick_TIFF,
                                                          rep_TIFF,
                                                          nb_channel), axis=1)

    print "IP extracted for Scan", cycle_number
    return IP_sub


def main_extractIP(read_noExtension_path, df_scan, path_save_IPextraction, transform_matrix_Rep,
                   transform_matrix_Nick, path_illumination, save=True):
    # Here the dataframe is grouped by Scan so that the TIFF is open once for the given Scan
    All_IP_Yoyo = df_scan.groupby(['Scan']).apply(IP_extraction,
                                                  transform_matrix_rep=transform_matrix_Rep,
                                                  transform_matrix_nick=transform_matrix_Nick,
                                                  path_illu=path_illumination,
                                                  path_no_extension=read_noExtension_path)
    print 'Intensity profile --> done'

    df_final_IPextraction = df_scan
    df_final_IPextraction['IP_Yoyo'] = All_IP_Yoyo.get(1)
    df_final_IPextraction['IP_Rep'] = All_IP_Yoyo.get(2)
    df_final_IPextraction['IP_Nick'] = All_IP_Yoyo.get(3)

    if save == True:
        df_final_IPextraction.to_csv(path_or_buf=path_save_IPextraction, sep='\t')

    return df_final_IPextraction, All_IP_Yoyo


df_allScans_path = "path/Molecules.mol"

print('\n Importing molecules: ' + df_allScans_path)
cycle_number = 1
df_mol = pd.read_csv(df_allScans_path, header=3, sep='\t')  # Importing molecules
df_mol.drop(df_mol.columns[len(df_mol.columns) - 1], axis=1, inplace=True)
df_Scan = df_mol[df_mol.Scan == cycle_number]
df_Scan['MoleculeSlope'] = (df_Scan['YEnd'] - df_Scan['YStart']) / (df_Scan['XEnd'] - df_Scan['XStart'])
print('--> Done.')

# Path
date = 'date'
scan_nbr = "scan_nbr"
run_name = "runname"
general_path = "path"
readimg_noExtension_path = "path"

transform_matrix = np.load('path')
transform_matrixN = np.load('path')
path_illu = 'path'

path_save_IPextraction = "/path" + date + ".csv"
result = main_extractIP(readimg_noExtension_path, df_Scan, path_save_IPextraction, transform_matrix,
                        transform_matrixN, path_illu, save=True)
