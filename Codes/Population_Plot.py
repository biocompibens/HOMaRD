from ast import literal_eval
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import misc


def plot_Population_Rep_avg(df_allmappedMol, df_rep, save_path, plot_img=False, threshold=False, save_file=True,
                            save_img=False, sorting='MoleculeStart_onRef'):
    df_sorted = df_rep.sort_values([sorting], ascending=[False])

    ip_to_plot = 'IP_Rep'
    population_image = np.zeros([len(df_allmappedMol), int(df_allmappedMol.RefLen.get_values()[0] / 553.5) / 2])
    coverage_image = np.zeros([len(df_allmappedMol), int(df_allmappedMol.RefLen.get_values()[0] / 553.5) / 2])
    pos_Y = 0
    mol_id = df_sorted.MoleculeId.get_values()
    print 'replication signal', len(mol_id), len(df_rep), len(df_sorted)

    for molecule_ in mol_id:
        line_molecule = df_sorted[df_sorted.MoleculeId == molecule_]
        if line_molecule.Orientation.get_values()[0] == '+':
            IP_molecule0 = line_molecule[ip_to_plot].get_values()[0]
            if threshold == True:
                IP_molecule = np.array(IP_molecule0) > 565
            else:
                IP_molecule = IP_molecule0

            for pos_X in range(0, len(IP_molecule)):
                population_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = IP_molecule[pos_X]
                coverage_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = 1
            pos_Y += 1
        else:
            IP_molecule0 = line_molecule[ip_to_plot].get_values()[0][::-1]
            if threshold == True:
                IP_molecule = np.array(IP_molecule0) > 565
            else:
                IP_molecule = IP_molecule0

            for pos_X in range(0, len(IP_molecule)):
                population_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = IP_molecule[pos_X]
                coverage_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = 1
            pos_Y += 1

    if plot_img == True:
        plt.figure()
        plt.imshow(population_image[:len(mol_id)][::-1])

    if save_img == True:
        scipy.misc.imsave(save_path + str(ip_to_plot) + '_PopulationPlot.tiff', population_image[:len(mol_id)])

    if save_file == True:
        np.save(save_path + 'Coverage_NbrSample' + str(len(mol_id)), coverage_image)
        np.save(save_path + 'Replication_NbrSample' + str(len(mol_id)), population_image[:len(mol_id)])

    print 'population image', np.shape(population_image[:len(mol_id)])
    print 'coverage image', np.shape(coverage_image[:len(mol_id)])
    print len(mol_id)
    return population_image[:len(mol_id)], mol_id, coverage_image[:len(mol_id)]


def plot_Population_Nick(df_allmappedMol, df_nick, mol_id, save_path, save_file=True, save_img=False, plot_img=False):
    # df_mapped_sup90_1FOV_0 = df_nick[df_nick.MoleculeStart_onRefINT >= 0]
    df_mapped_sup90_1FOV_0 = df_nick
    print len(mol_id), len(df_nick)

    ip_to_plot = 'IP_Nick'
    population_image = np.zeros([len(df_allmappedMol), int(df_allmappedMol.RefLen.get_values()[0] / 553.5) / 2])
    population_image_yoyo = np.zeros([len(df_allmappedMol), int(df_allmappedMol.RefLen.get_values()[0] / 553.5) / 2])
    pos_Y = 0

    for molecule_ in mol_id:
        line_molecule = df_mapped_sup90_1FOV_0[df_mapped_sup90_1FOV_0.MoleculeId == molecule_]
        if line_molecule.Orientation.get_values()[0] == '+':
            IP_molecule = line_molecule[ip_to_plot].get_values()[0]
            IP_molecule_Yoyo = line_molecule['IP_Yoyo'].get_values()[0]
            for pos_X in range(0, len(IP_molecule)):
                population_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = IP_molecule[pos_X]
                population_image_yoyo[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = \
                IP_molecule_Yoyo[pos_X]
            pos_Y += 1
        else:
            IP_molecule = line_molecule[ip_to_plot].get_values()[0][::-1]
            IP_molecule_Yoyo = line_molecule['IP_Yoyo'].get_values()[0][::-1]
            for pos_X in range(0, len(IP_molecule)):
                population_image[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = IP_molecule[pos_X]
                population_image_yoyo[pos_Y][pos_X + int(np.round(line_molecule.MoleculeStart_onRef))] = \
                IP_molecule_Yoyo[pos_X]
            pos_Y += 1

    if plot_img == True:
        plt.figure()
        plt.imshow(population_image[:len(mol_id)][::-1])

    if save_img == True:
        scipy.misc.imsave(save_path + str(ip_to_plot) + '_PopulationPlot.tiff', population_image[:len(mol_id)])
        scipy.misc.imsave(save_path + str('IP_Yoyo') + '_PopulationPlot.tiff', population_image_yoyo[:len(mol_id)])

    if save_file == True:
        np.save(save_path + 'Nick_NbrSample' + str(len(mol_id)), population_image[:len(mol_id)])
        np.save(save_path + 'Yoyo_NbrSample' + str(len(mol_id)), population_image_yoyo[:len(mol_id)])

    return population_image[:len(mol_id)], population_image_yoyo[:len(mol_id)]


def plot_populationImage(replication_Pop_img, nick_Pop_img, yoyo_Pop_img, title, path_save, save_img=True):
    fig2 = plt.figure(figsize=(12, 15))
    plt.suptitle(title)
    plt.subplot(131)
    plt.imshow(replication_Pop_img[::-1])
    plt.subplot(132)
    plt.imshow(nick_Pop_img[::-1])
    plt.subplot(133)
    plt.imshow(yoyo_Pop_img[::-1])
    if save_img == True:
        fig2.savefig(path_save + 'PopulationImage.svg', dpi=500)


def main_population_plot(replication_Pop_img, nick_Pop_img, coverage_img, plottitle_, path_save, save_img=True,
                         save_files=True):
    coverage_count = np.count_nonzero(coverage_img, axis=0)

    # IP Rep
    sum_rep = np.sum(replication_Pop_img, axis=0)
    norm_sum_rep = np.divide(sum_rep, coverage_count)
    figure_rep = plot_population_profile(coverage_count, sum_rep, norm_sum_rep, plottitle_, color='DeepPink',
                                         channel='Replication')

    # IP Nick
    sum_nick = np.sum(nick_Pop_img, axis=0)
    norm_sum_nick = np.divide(sum_nick, coverage_count)
    figure_nick = plot_population_profile(coverage_count, sum_nick, norm_sum_nick, plottitle_, color='LimeGreen',
                                          channel='Nick')

    if save_img == True:
        figure_rep.savefig(path_save + 'Coverage' + str(coverage_img.shape[0]) + 'Replication.png', dpi=250)
        figure_nick.savefig(path_save + 'Coverage' + str(coverage_img.shape[0]) + 'Nick.png', dpi=250)

    if save_files == True:
        np.save(path_save + 'coverage', coverage_count)
        np.save(path_save + 'sum_rep', sum_rep)
        np.save(path_save + 'norm_sum_rep', norm_sum_rep)
        np.save(path_save + 'sum_nick', sum_nick)
        np.save(path_save + 'norm_sum_nick', norm_sum_nick)

    return coverage_count, sum_rep, norm_sum_rep, sum_nick, norm_sum_nick


def plot_population_profile(coverage, sum, normsum, title, color='DeepPink', channel='Replication'):
    fig1 = plt.figure(figsize=(20, 15))
    plt.suptitle(title)
    plt.subplot(311)
    plt.title('Coverage')
    plt.plot(coverage, color='Gray')
    plt.xlim([0, 1200])

    plt.subplot(312)
    plt.title('Sum ' + channel + ' IP')
    plt.plot(sum, color=color)
    plt.xlim([0, 1200])

    plt.subplot(313)
    plt.title('Sum/Coverage ' + channel + ' IP')
    plt.plot(normsum, color=color)
    plt.xlim([0, 1200])
    return fig1


date = 'date'
input_IP_root = "path"
mol_allmapped_Corr_path = input_IP_root + 'path'  # All mapped molecules (1FOV + otherFOV)
mol_mapped_Corr_path = 'path'

df_mapped_sup90_1FOV_replicated_mean = pd.read_csv(mol_mapped_Corr_path,
                                                   sep='\t')  # Molecules that have mapped and are in pop HIGH
Mol_ID = df_mapped_sup90_1FOV_replicated_mean.MoleculeId.get_values()
df_selected = pd.DataFrame({'MoleculeId': Mol_ID})

df_mappedMol_all = pd.read_csv(mol_allmapped_Corr_path, sep='\t')
df_mappedMol_all.drop(df_mappedMol_all.columns[0], axis=1, inplace=True)
df_mappedMol = df_mappedMol_all[df_mappedMol_all.IsStitched == 0]

df_mapped_sup90_1FOV_replicated_mean_contig1 = pd.merge(df_selected, df_mappedMol, on=['MoleculeId'])

df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Rep = df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Rep.apply(
    literal_eval)
df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Nick = df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Nick.apply(
    literal_eval)
df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Yoyo = df_mapped_sup90_1FOV_replicated_mean_contig1.IP_Yoyo.apply(
    literal_eval)

df_mapped_sup90_1FOV_replicated_mean_contig1_MolStartOnrefSUP0 = df_mapped_sup90_1FOV_replicated_mean_contig1[
    df_mapped_sup90_1FOV_replicated_mean_contig1.MoleculeStart_onRef > 0]

save_path_file_all8 = "path" + date + "/"

img_rep_allmol_mean = plot_Population_Rep_avg(df_mappedMol,
                                              df_mapped_sup90_1FOV_replicated_mean_contig1_MolStartOnrefSUP0,
                                              save_path_file_all8, save_img=True, plot_img=False, save_file=True)
img_nick_allmol_mean = plot_Population_Nick(df_mappedMol,
                                            df_mapped_sup90_1FOV_replicated_mean_contig1_MolStartOnrefSUP0,
                                            img_rep_allmol_mean[1],
                                            save_path_file_all8, save_img=True, plot_img=False, save_file=True)

title_ = 'Replicated molecules in Contig1' + str(img_rep_allmol_mean[0].shape[0])
plot_img = plot_populationImage(img_rep_allmol_mean[0], img_nick_allmol_mean[0], img_nick_allmol_mean[1], title_,
                                save_path_file_all8, save_img=True)
a = main_population_plot(img_rep_allmol_mean[0], img_nick_allmol_mean[0], img_rep_allmol_mean[2],
                         'MolStartonRef>0 1FOV', save_path_file_all8, save_img=True)
