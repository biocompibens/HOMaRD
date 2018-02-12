from ast import literal_eval
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random


def split_channels_Illu(img_Tiff, cycle_number, path_illu):  # Read and split the TIFF image
    import tifffile
    yoyo_channel = []
    rep_channel = []
    nick_channel = []
    if cycle_number < 10:
        img_path_init = (img_Tiff + str("00") + str(cycle_number) + ".tiff")
    else:
        img_path_init = (img_Tiff + str("0") + str(cycle_number) + ".tiff")

    original_images_1140 = tifffile.imread(img_path_init)  # read the TIFF image
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


#
# ############################################### INPUT/OUTPUT paths ###################################################

xmap_merge = "path"
nick_pos_trans_path = 'path'
transform_matrix_Replication = np.load('path')
transform_matrix_Nick = np.load('path')
illumination_fct_path = 'path'

date = 'date'
run_name = "run_name"
input_mapping_root = "path"
ref_NickPos_path = input_mapping_root + "MoleculeQualityReport_r.cmap"
scan_rescaleFact_path = input_mapping_root + "MoleculeQualityReport.scan"
global_rescaleFact_path = input_mapping_root + "MoleculeQualityReport.txt"

############ LOADING FILES ############

Contig_val = 1
df_mappedMol_ori = pd.read_csv(xmap_merge, sep='\t')  # Merged IP + Xmap mapped molecules
df_mappedMol0 = df_mappedMol_ori[df_mappedMol_ori.IsStitched == 0]  # Molecule only on 1 field of view
df_mappedMol0 = df_mappedMol0[
    df_mappedMol0.RefContigID == Contig_val]  # choose only molecules that are on a given contig/chromosome

df_labs = pd.read_csv(nick_pos_trans_path, sep='\t')
df_ref_nick = pd.read_csv(ref_NickPos_path, header=8, skiprows=[9, 9], sep='\t')  # Nick position on reference genome
df_ref_nick_contig = df_ref_nick[df_ref_nick['#h CMapId'] == Contig_val]
df_rescale_scan = pd.read_csv(scan_rescaleFact_path, sep='\t', header=4)  # Scaling factor (scan depend)
general_SF = []  # Chip dependant rescaling factor

for line in open(global_rescaleFact_path, 'r'):
    if "Bpp" in line:
        general_SF.append(float(line.split()[1]) / 500)

Scan_number = 5
df_mappedMol = df_mappedMol0.groupby(df_mappedMol0.Scan).get_group(
    Scan_number)  # Molecules that are contained in the given Scan
molID = df_mappedMol['MoleculeId'].get_value(
    random.SystemRandom().choice(df_mappedMol.index))  # Select random ID in the given dataframe

# MOLECULE SELECTION
line_molecule = df_mappedMol[
    df_mappedMol['MoleculeId'] == molID]  # Get the line corresponding to the molecule ID ALL CORR
scan_nb = line_molecule['Scan'].get_values()[0]  # Get scan number
fov_nb = line_molecule['FovStart'].get_values()[0]  # Get FOV

# SCALING FACTOR
scan_scalingfact = df_rescale_scan[df_rescale_scan['Scan'] == (scan_nb - 1)].get_values()[0][3]
finalScalingFactor = general_SF[0] * scan_scalingfact

# MOLECULE and REFERENCE LABELS
labels_positions0 = df_labs[
    df_labs.MoleculeID == molID]  # Get the labels from molecules in the .lab file that has the transformed positions
labels_positions = labels_positions0[
    labels_positions0.SNR > 2.75]  # Filter labels on SNR value according to the filter used for the mapping
allRef_positions = (df_ref_nick_contig.Position.get_values() / (float(
    553.5) * finalScalingFactor))  # All positions of the reference with scaling factor applied and converted in basepair

# INTENSITY PROFILES with and without correction
molecule_orientation = line_molecule.Orientation.get_values()[0]
ip_mol_yoyo = literal_eval(line_molecule.IP_Yoyo.get_values()[0])
ip_mol_rep = literal_eval(line_molecule.IP_Rep.get_values()[0])
ip_mol_nick = literal_eval(line_molecule.IP_Nick.get_values()[0])

########## CHECK LABEL POSITION ON EXTRACTED IP (Nick) ##############

# Molecule IP and transformed nick positions plot : we check here if the transformed labels are well positionned on the
# extracted intensity profile
save_check_nickPos = False
transformed_positions = labels_positions.PositionOnMolTrans.get_values()  # Transformed position of the labels found on
# molecules
fig_check_position = plt.figure(figsize=(20, 5))
plt.plot(ip_mol_nick, color='limeGreen')  # Nick IP with both corrections
for pos_on_mol in transformed_positions:
    plt.axvline(pos_on_mol, color='green')
# Plotting the position of the first well mapped nick
first_nick_mapped_MOL0 = line_molecule.Qry_firstNick.get_values()[0]
plt.axvline(transformed_positions[first_nick_mapped_MOL0 - 1], color='r')
if save_check_nickPos == True:
    fig_check_position.savefig('fig_path' + str(Scan_number) + 'molID' + str(molID) + '_LabelDetection.png')

########## MAPPED MOLECULE VISU ##############

# Creating other variables having the IP of each color
NickIP_onreference = ip_mol_nick
RepIP_onreference = ip_mol_rep
YoyoIP_onreference = ip_mol_yoyo

MoleculeStart_onRef = line_molecule.MoleculeStart_onRef.get_values()[0]
x = np.arange(0, len(
    NickIP_onreference)) + MoleculeStart_onRef  # Position the IP on the reference genome from the 1st well mapped nick

if molecule_orientation == '+':
    yN, yR, yY = NickIP_onreference, RepIP_onreference, YoyoIP_onreference
else:
    yN, yR, yY = NickIP_onreference[::-1], RepIP_onreference[::-1], YoyoIP_onreference[::-1]

length_noScaling = line_molecule['Length'].get_values()[
    0]  # Length of the molecule in the .mol file (Autodetect output)
length_extracted = len(literal_eval(line_molecule['IP_Yoyo'].get_values()[0]))
length_expected_Scaling = line_molecule['Length'].get_values()[0] * finalScalingFactor
length_xmap = line_molecule['QryLen'].get_values()[0]  # Length of the mapped molecule
cutoff = line_molecule['IsCutOff'].get_values()[0]
df_length = pd.DataFrame({"1.Original molecule (.mol) ": [length_noScaling], "2.IP (px) ": [length_extracted],
                          "4.Mapped molecule (px) ": [length_xmap / 553.5],
                          "3.Expected (px)": [length_expected_Scaling]},
                         index=['Length'])  # after resampling (scaling factor applied)

fig1 = plt.figure(figsize=(20, 15))
plt.suptitle('Single molecule intensity profile (in px)' + '\n' + ' Molecule ID: ' + str(
    line_molecule['MoleculeId'].get_values()[0]) + ','
             + ' Confidence: ' + str(line_molecule['Confidence'].get_values()[0]) + ','
             + ' Orientation: ' + str(line_molecule['Orientation'].get_values()[0]) + ','
             + ' CutOff: ' + str(cutoff) + ','
             + ' AutodetectLen :' + str(np.round(length_noScaling)) + ','
             + ' ExtractedLen :' + str(length_extracted) + ','
             + ' MappedLen :' + str(np.round(length_xmap) / float(553.5)) + ','
             + '\n'
             + ' ResampledLenAll :' + str(len(ip_mol_nick) * finalScalingFactor),
             fontsize=16)

ax = plt.subplot(311)
plt.plot(x, yY, color='DodgerBlue', label='Corrected IP', alpha=1, lw=2)
plt.xlim([x[0] - 5, x[-1] + 5])
plt.legend()

ax3 = plt.subplot(312)
plt.plot(x, yR, color='DeepPink', label='Corrected IP', alpha=1, lw=2)
plt.xlim([x[0] - 5, x[-1] + 5])
plt.legend()

ax5 = plt.subplot(313)
# all reference positions
for ref_allpos in allRef_positions:
    plt.axvline(ref_allpos, color='gray', alpha=0.5)
ref_IDpos_lineMol = literal_eval(line_molecule.Ref_AllNicks.get_values()[0])
for ref_pos in ref_IDpos_lineMol:
    ref_pos_lineMol = ((df_ref_nick_contig.Position.get_value(ref_pos - 1) / (
    float(553.5) * finalScalingFactor)))
    plt.axvline(ref_pos_lineMol, color='black', alpha=0.5)

# Plot of the IP
plt.plot(x, yN, color='LimeGreen', label='Corrected IP', alpha=1, lw=2)
plt.xlim([x[0] - 5, x[-1] + 5])
plt.legend()
