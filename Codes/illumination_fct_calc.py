import numpy as np
import tifffile


def split_flip_img(path_NOextension, cycle_number):  # Read and split the TIFF image

    yoyo_channel = []
    rep_channel = []
    nick_channel = []

    if cycle_number < 10:
        img_path_init = (path_NOextension + '_Scan' + str("00") + str(cycle_number) + ".tiff")
    else:
        img_path_init = (path_NOextension + '_Scan' + str("0") + str(cycle_number) + ".tiff")

    original_images_1140 = tifffile.imread(img_path_init)  # read the TIFF image
    nb_channels = (len(original_images_1140) / 1140)  # Count the number of channels in the TIFF image
    # Flip the original image horizontally and vertically
    if nb_channels == 2:  # If the number of channel is 2
        for i in range(0, len(original_images_1140), 2):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])))
            rep_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])))
        print "Splitting Done --> 2 channels", " Scan", cycle_number
        return nb_channels, yoyo_channel, nick_channel

    if nb_channels == 3:  # If the number of channel is 3
        for i in range(0, len(original_images_1140), 3):
            yoyo_channel.append(np.flipud(np.fliplr(original_images_1140[i])))
            rep_channel.append(np.flipud(np.fliplr(original_images_1140[i + 1])))
            nick_channel.append(np.flipud(np.fliplr(original_images_1140[i + 2])))
        print "Splitting Done --> 3 channels ", " Scan", cycle_number
        return nb_channels, yoyo_channel, rep_channel, nick_channel


def med_calc_run(path_img_NOextension):
    med_img_per_cycleYOYO = []
    med_img_per_cycleREP = []
    med_img_per_cycleNICK = []

    for cycle_number in range(1, 31):
        splitted_image = split_flip_img(path_img_NOextension, cycle_number=cycle_number)
        print "--> Image Splitted ", cycle_number

        if splitted_image[0] == 3:
            yoyoTiff = splitted_image[1]
            repTiff = splitted_image[2]
            nickTiff = splitted_image[3]
            yoyo_med_img = np.median(yoyoTiff, axis=0)
            rep_med_img = np.median(repTiff, axis=0)
            nick_med_img = np.median(nickTiff, axis=0)
            med_img_per_cycleYOYO.append(yoyo_med_img)
            med_img_per_cycleREP.append(rep_med_img)
            med_img_per_cycleNICK.append(nick_med_img)
        else:
            yoyoTiff = splitted_image[1]
            nickTiff = splitted_image[2]
            yoyo_med_img = np.median(yoyoTiff, axis=0)
            nick_med_img = np.median(nickTiff, axis=0)
            med_img_per_cycleYOYO.append(yoyo_med_img)
            med_img_per_cycleREP.append(0)
            med_img_per_cycleNICK.append(nick_med_img)

        print 'Scan', cycle_number, '--> Done'
    return med_img_per_cycleYOYO, med_img_per_cycleREP, med_img_per_cycleNICK


general_path = "path"
save_path = 'path'
run_name = "run_name"
date = 'date'
path_img_NOextension = general_path + run_name
all_colors_med_image_cycle = med_calc_run(path_img_NOextension)
save_file = np.save(save_path + 'IlluminationFct_' + run_name + '_' + date, all_colors_med_image_cycle)
