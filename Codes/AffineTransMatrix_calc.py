import matplotlib
matplotlib.use('Agg')
import numpy as np
import tifffile
import scipy
from scipy import optimize
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
import skimage
from skimage import morphology
from scipy.ndimage import filters


def split_flip_img(img_Tiff):
    # Read and split the TIFF image
    yoyo_channel = []
    nick_channel = []
    rep_channel = []

    original_images_1140 = tifffile.imread(img_Tiff)  # read the TIFF image
    nb_channels = (len(original_images_1140) / 1140)  # Count the number of channels in the TIFF image

    # Flip the original image horizontally and vertically
    if nb_channels == 2:  # If the number of channel is 2
        for i in range(0, len(original_images_1140), 2):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])))
            nick_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])))
        print "FLIP Splitting Done 2 channels"
        return 2, yoyo_channel, nick_channel
    if nb_channels == 3:  # If the number of channel is 3
        for i in range(0, len(original_images_1140), 3):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])))
            rep_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])))
            nick_channel.append(np.flipud(np.fliplr(original_images_1140[i + 2])))
        print "FLIP Splitting Done 3 channels "
        return 3, yoyo_channel, rep_channel, nick_channel


def Local_max_REP(original_image, rect_erosion, rect_locmax, threshold_value, kernel_channel):
    # Creating mask
    erosion = skimage.morphology.erosion(original_image, selem=rect_erosion)
    erosion_subtract = original_image - erosion
    diff_thresh = (erosion_subtract > np.percentile(erosion_subtract, threshold_value))
    dilation = skimage.morphology.opening(diff_thresh, selem=kernel_channel)
    masked_img = original_image * dilation

    # Keeping only the maximum contained in the mask
    final_molRep_max = skimage.morphology.local_maxima(masked_img, selem=rect_locmax)

    return erosion, erosion_subtract, diff_thresh, dilation, masked_img, final_molRep_max


def Local_max_Nick(original_image, rect_erosion, threshold_value, kernel_channel, neighborhood_size_value):
    # Detecting local max on original image
    Max_originalImg = filters.maximum_filter(original_image, neighborhood_size_value)
    All_max_ori = (original_image == Max_originalImg)
    # Creating mask
    erosion = skimage.morphology.erosion(original_image, selem=rect_erosion)
    erosion_subtract = original_image - erosion
    diff_thresh = (erosion_subtract > np.percentile(erosion_subtract, threshold_value))
    dilation = skimage.morphology.opening(diff_thresh, selem=kernel_channel)
    masked_img = original_image * dilation
    # Keeping only the maximum contained in the mask
    All_max_ori[masked_img == 0] = 0
    final_molNick_max = All_max_ori
    return erosion, erosion_subtract, diff_thresh, dilation, masked_img, final_molNick_max


def Local_max_Yoyo(original_image, rect_erosion, rect_locmax, threshold_value, kernel_channel):
    # Creating mask
    erosion = skimage.morphology.erosion(original_image, selem=rect_erosion)
    erosion_subtract = original_image - erosion
    diff_thresh = (erosion_subtract > np.percentile(erosion_subtract, threshold_value))
    dilation = skimage.morphology.dilation(diff_thresh, selem=kernel_channel)
    # Detecting local max on original image
    masked_img = original_image * dilation
    final_molYoyo_max = skimage.morphology.local_maxima(masked_img, selem=rect_locmax)
    return erosion, erosion_subtract, diff_thresh, dilation, masked_img, final_molYoyo_max


def registration_REP(local_max_ref, local_max_work, ref_image_frame, verbose=False):
    max_yoyo0 = np.swapaxes(np.argwhere(local_max_ref[5] == 1), 1, 0)
    max_rep0 = np.swapaxes(np.argwhere(local_max_work[5] == 1), 1, 0)

    locmax_yoyo = []
    locmax_rep = []
    for indloc in range(0, len(max_yoyo0[0])):
        locmax_yoyo.append([max_yoyo0[1][indloc], max_yoyo0[0][indloc]])
    for j in range(0, len(max_rep0[0])):
        locmax_rep.append([max_rep0[1][j], max_rep0[0][j]])

    dist_matrix = dist.cdist(locmax_rep, locmax_yoyo)
    index_min_dist = np.argmin(dist_matrix, axis=1)

    data_ref = []
    data_work = []
    k = 0
    for ind in index_min_dist:
        if dist_matrix[k][ind] < 6:
            data_ref.append([locmax_yoyo[ind][0], locmax_yoyo[ind][1], 1])
            data_work.append([locmax_rep[k][0], locmax_rep[k][1], 1])
        k += 1
    ref_ = np.array(data_ref)
    work_ = np.array(data_work)

    if verbose:
        fig3 = plt.figure()
        for i in range(0, len(ref_)):
            plt.plot((ref_[i][0], work_[i][0]), (ref_[i][1], work_[i][1]))
        plt.imshow(ref_image_frame, cmap='gray')

    init_matrix = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1])
    other_value = 5
    scal_val = 0.99
    parameter_bounds = [[0.99, 0.99], [-other_value, other_value], [0, 0],
                        [-other_value, other_value], [0.99, 0.99], [0, 0],
                        [-other_value, other_value], [-other_value, other_value], [1, 1]]
    transform_mat_found = scipy.optimize.minimize(transform_func, init_matrix, args=(ref_, work_, scal_val),
                                                  bounds=parameter_bounds, method='SLSQP')

    return transform_mat_found


def registration_NICK(local_max_ref, local_max_work, ref_image_frame, verbose=False):
    max_yoyo0 = np.swapaxes(np.argwhere(local_max_ref[5] == 1), 1, 0)
    max_rep0 = np.swapaxes(np.argwhere(local_max_work[5] == 1), 1, 0)

    locmax_yoyo = []
    locmax_rep = []
    for indloc in range(0, len(max_yoyo0[0])):
        locmax_yoyo.append([max_yoyo0[1][indloc], max_yoyo0[0][indloc]])
    for j in range(0, len(max_rep0[0])):
        locmax_rep.append([max_rep0[1][j], max_rep0[0][j]])

    dist_matrix = dist.cdist(locmax_rep, locmax_yoyo)
    index_min_dist = np.argmin(dist_matrix, axis=1)

    data_ref = []
    data_work = []
    k = 0
    for ind in index_min_dist:
        if dist_matrix[k][ind] < 6:
            data_ref.append([locmax_yoyo[ind][0], locmax_yoyo[ind][1], 1])
            data_work.append([locmax_rep[k][0], locmax_rep[k][1], 1])
        k += 1
    ref_ = np.array(data_ref)
    work_ = np.array(data_work)

    if verbose:
        fig3 = plt.figure()
        for i in range(0, len(ref_)):
            plt.plot((ref_[i][0], work_[i][0]), (ref_[i][1], work_[i][1]))
        plt.imshow(ref_image_frame, cmap='gray')

    init_matrix = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1])
    other_value = 5
    scal_val = 0.996
    parameter_bounds = [[0.996, 0.996], [-other_value, other_value], [0, 0],
                        [-other_value, other_value], [0.996, 0.996], [0, 0],
                        [-other_value, other_value], [-other_value, other_value], [1, 1]]

    transform_mat_found = scipy.optimize.minimize(transform_func, init_matrix, args=(ref_, work_, scal_val),
                                                  bounds=parameter_bounds, method='SLSQP')
    return transform_mat_found


def transform_func(transformation_matrix, ref_data, work_data, scaling_val):
    p1, p2, p3, p4, p5, p6, p7, p8, p9 = transformation_matrix
    p1 = scaling_val
    p2 = 0
    p4 = 0
    matrix = np.array([[p1, p2, p3],
                       [p4, p1, p6],
                       [p7, p8, p9]])
    transformed_coord = np.dot(work_data, matrix)
    cost_val = cost_calc(ref_data, transformed_coord)
    return cost_val


def cost_calc(ref_data, transformed_coord_work):
    k = 0
    cost = 0
    dist_matrix_trans = dist.cdist(ref_data, transformed_coord_work)
    index_min_dist_trans = np.argmin(dist_matrix_trans, axis=1)
    for i in index_min_dist_trans:
        cost += dist_matrix_trans[k][i]
        k += 1
    return cost


def register_all_image(image_path, proportion):
    num, yoyo_, rep_, nick_ = split_flip_img(image_path)  # split the TIFF image

    # Normalisation of the images
    all_avg_YOYO = []
    for j in range(0, 1140):
        all_avg_YOYO.append(np.mean(yoyo_[j]))
    all_avg_sort = np.argsort(
        all_avg_YOYO)  # sort the average_YOYO and return the position for each img in the same order
    selected_images = all_avg_sort[(len(all_avg_YOYO) * proportion / 100):]  # only keep 20% of the images (once sorted)
    print "Selected Images : ", len(selected_images), "/1140"

    AVG_image_yoyo = []
    transform_per_image_rep = []
    transform_per_image_nick = []
    count = 0
    for i in selected_images:
        # select the image in "sel_img array
        ref_img = yoyo_[i]
        work_img_rep = rep_[i]
        work_img_nick = nick_[i]

        # calculate the average value of the image YOYO for the selected images
        AVG_image_yoyo.append(np.average(ref_img))

        # Common kernels
        rect_erosion = np.array([[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]], dtype='uint8')
        rect_locmax = np.array([[0, 0, 0, 0, 0], [1, 1, 1, 1, 1], [0, 0, 0, 0, 0]], dtype='uint8')
        # Channel specific kernel and value
        disk_opening = skimage.morphology.selem.disk(1)
        rect_dilation_YOYO = skimage.morphology.selem.rectangle(5, 2)
        threshold_value = 90
        threshold_value_YOYO = 85

        local_max_ref = Local_max_Yoyo(ref_img, rect_erosion, rect_locmax, threshold_value=threshold_value_YOYO,
                                       kernel_channel=rect_dilation_YOYO)
        local_max_work_rep = Local_max_REP(work_img_rep, rect_erosion, rect_locmax, threshold_value=threshold_value,
                                           kernel_channel=disk_opening)
        local_max_work_nick = Local_max_Nick(work_img_nick, rect_erosion, threshold_value=threshold_value,
                                             kernel_channel=disk_opening, neighborhood_size_value=5)

        transform_matrix_rep = registration_REP(local_max_ref, local_max_work_rep, ref_img)
        transform_matrix_nick = registration_NICK(local_max_ref, local_max_work_nick, ref_img)
        transform_per_image_rep.append(transform_matrix_rep.x)
        transform_per_image_nick.append(transform_matrix_nick.x)
        count += 1

    print "DONE"
    return transform_per_image_rep, transform_per_image_nick, selected_images, AVG_image_yoyo


def main(path_NOextension, scanStart, scanEndIncluded, verbose=False):
    def med_matrix(final_matrices):
        p1, p2, p3, p4, p5, p6, p7, p8, p9 = [], [], [], [], [], [], [], [], []
        x_axis = []
        all_matrices = []

        for i in range(0, len(final_matrices)):
            x_axis.append(i)
            p1.append(final_matrices[i][0])
            p2.append(final_matrices[i][1])
            p3.append(final_matrices[i][2])
            p4.append(final_matrices[i][3])
            p5.append(final_matrices[i][4])
            p6.append(final_matrices[i][5])
            p7.append(final_matrices[i][6])
            p8.append(final_matrices[i][7])
            p9.append(final_matrices[i][8])
        all_matrices.append([np.median(p1), np.median(p2), np.median(p3),
                             np.median(p4), np.median(p5), np.median(p6),
                             np.median(p7), np.median(p8), np.median(p9)])
        return all_matrices

    for Scan in range(scanStart, scanEndIncluded + 1):

        if Scan < 10:
            img_path_init = (path_NOextension + '_Scan' + str("00") + str(Scan) + ".tiff")
            print img_path_init
        else:
            img_path_init = (path_NOextension + '_Scan' + str("0") + str(Scan) + ".tiff")
            print img_path_init

        final_matrices_rep, final_matrices_nick, image_number, average_value_images = register_all_image(img_path_init,
                                                                                                         proportion=85)

        replication_final_mat = med_matrix(final_matrices_rep)
        nick_final_mat = med_matrix(final_matrices_nick)

    return final_matrices_rep, final_matrices_nick, image_number, replication_final_mat, nick_final_mat


save_file = True
run_name = "path"
path_noext = "path"
save_path = 'path'
date = 'date'

all_mat = main(path_noext, 1, 2)
if save_file == True:
    np.save(save_path + 'AffineTransMat_RED_' + run_name + '_' + date, all_mat[0])
    np.save(save_path + 'AffineTransMat_NICK_' + run_name + '_' + date, all_mat[1])
    np.save(save_path + 'AffineTransMat_selectedImg_' + run_name + '_' + date, all_mat[2])
    np.save(save_path + 'AffineTransMat_RED_median_' + run_name + '_' + date, all_mat[3])
    np.save(save_path + 'AffineTransMat_NICK_median_' + run_name + '_' + date, all_mat[4])
