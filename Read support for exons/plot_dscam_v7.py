import glob
import os
import pandas as pd
import svgwrite
from svgwrite import mm

'''
Description: require python >= 3.6
height bar is read counts (normalized); blue is exon; green is known donor; yellow is known acceptor; purple is novel donor; orange is novel acceptor.
'''

# ==========parameters to be changed==========
# ==========parameters for adjusting images==========
xoffset = 30
yoffset = 30
exon_height = 3
xaxis_length = 150
yaxis_length = 150
# ==========parameters for adjusting images==========
min_junction_reads_threshold = 10
splice_junction_file_suffix = 'SJ.out.tab'
wanted = 'FBgn0033159'  # wanted string which is contained in last column of the GTF file
path_gtf_file = '/home/user/data/dscam_fruit_fly.gtf'  # can be a plain text or compressed file (.gz, .zip, .xz, .bz2)
path_splice_junction_folder = '/home/user/data/fruit_fly_sj_data'
path_output_svg_folder = '/home/user/data/fruit_fly_plotting'
# ==========parameters to be changed==========


def draw_sj_image(df_gtf, df_splice_junction, path_output_svg_file, img_title):
    df_exon_drawing = df_gtf[['start', 'end', 'info']].copy()
    df_exon_drawing['exon_number'] = df_exon_drawing['info'].str.extract(r'exon_number "(.*?)";')
    df_exon_drawing['exon_number_count'] = df_exon_drawing['exon_number'].map(df_exon_drawing['exon_number'].value_counts())
    df_exon_drawing['exon_number_cumulative_count'] = df_exon_drawing.groupby('exon_number').cumcount() + 1    # need to add 1
    df_exon_drawing['exon_numering'] = df_exon_drawing['exon_number']
    df_exon_drawing.loc[df_exon_drawing['exon_number_count'] > 1, 'exon_numering'] = df_exon_drawing.loc[df_exon_drawing['exon_number_count'] > 1, 'exon_number'].astype('str') + '.' + df_exon_drawing.loc[df_exon_drawing['exon_number_count'] > 1, 'exon_number_cumulative_count'].astype('str')

    min_coordinate = df_exon_drawing['start'].min()
    max_coordinate = df_exon_drawing['end'].max()
    global xoffset, yoffset, exon_height, xaxis_length, yaxis_length
    df_exon_drawing['scaled_start'] = (df_exon_drawing['start'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_exon_drawing['scaled_end'] = (df_exon_drawing['end'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_exon_drawing['scaled_width'] = df_exon_drawing['scaled_end'] - df_exon_drawing['scaled_start']
    # draw exons
    dwg = svgwrite.Drawing(path_output_svg_file, profile='tiny')
    for idx_label in df_exon_drawing.index:
        dwg.add(dwg.text(df_exon_drawing.loc[idx_label, 'exon_numering'], insert=(df_exon_drawing.loc[idx_label, 'scaled_start'] + xoffset, yoffset - 2), fill='black', font_size="3.5px"))
        dwg.add(dwg.rect(insert=(df_exon_drawing.loc[idx_label, 'scaled_start'] + xoffset, yoffset), size=(df_exon_drawing.loc[idx_label, 'scaled_width'], exon_height), stroke='black', stroke_width=0.1, fill='lightblue', opacity=0.7))
    # process data for drawing splice junction reads
    df_sj_drawing = df_splice_junction[['start', 'end', 'number of uniquely mapping reads crossing the junction']].copy()
    df_sj_drawing['donor_exon_end'] = df_sj_drawing['start'] - 1
    df_sj_drawing['acceptor_exon_start'] = df_sj_drawing['end'] + 1
    df_sj_drawing.loc[df_sj_drawing['donor_exon_end'].isin(df_exon_drawing['end']), 'is_known_donor'] = True
    df_sj_drawing['is_known_donor'] = df_sj_drawing['is_known_donor'].fillna(False)
    df_sj_drawing.loc[df_sj_drawing['acceptor_exon_start'].isin(df_exon_drawing['start']), 'is_known_acceptor'] = True
    df_sj_drawing['is_known_acceptor'] = df_sj_drawing['is_known_acceptor'].fillna(False)
    df_temp_donor = pd.merge(df_exon_drawing[['start', 'end']], df_sj_drawing[['donor_exon_end']].reset_index(), how='right', left_on=['end'], right_on=['donor_exon_end'])
    df_temp_donor = df_temp_donor.rename({'start': 'donor_exon_start'}, axis=1)
    df_temp_donor = df_temp_donor.set_index('index').drop(['end', 'donor_exon_end'], axis=1)
    df_sj_drawing = pd.concat([df_sj_drawing, df_temp_donor], axis=1)
    df_temp_acceptor = pd.merge(df_exon_drawing[['start', 'end']], df_sj_drawing[['acceptor_exon_start']].reset_index(), how='right', left_on=['start'], right_on=['acceptor_exon_start'])
    df_temp_acceptor = df_temp_acceptor.rename({'end': 'acceptor_exon_end'}, axis=1)
    df_temp_acceptor = df_temp_acceptor.set_index('index').drop(['start', 'acceptor_exon_start'], axis=1)
    df_sj_drawing = pd.concat([df_sj_drawing, df_temp_acceptor], axis=1)

    df_sj_drawing['scaled_donor_exon_start'] = (df_sj_drawing['donor_exon_start'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_sj_drawing['scaled_donor_exon_end'] = (df_sj_drawing['donor_exon_end'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_sj_drawing['scaled_donor_exon_width'] = df_sj_drawing['scaled_donor_exon_end'] - df_sj_drawing['scaled_donor_exon_start']
    df_sj_drawing['scaled_acceptor_exon_start'] = (df_sj_drawing['acceptor_exon_start'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_sj_drawing['scaled_acceptor_exon_end'] = (df_sj_drawing['acceptor_exon_end'] - min_coordinate) / (max_coordinate - min_coordinate) * xaxis_length
    df_sj_drawing['scaled_acceptor_exon_width'] = df_sj_drawing['scaled_acceptor_exon_end'] - df_sj_drawing['scaled_acceptor_exon_start']

    df_sj_drawing['donor_exon_reads'] = df_sj_drawing.groupby('donor_exon_end')[['number of uniquely mapping reads crossing the junction']].transform('sum')
    df_sj_drawing['acceptor_exon_reads'] = df_sj_drawing.groupby('acceptor_exon_start')[['number of uniquely mapping reads crossing the junction']].transform('sum')

    total_reads = df_splice_junction['number of uniquely mapping reads crossing the junction'].sum()
    df_sj_drawing['scaled_donor_exon_reads'] = df_sj_drawing['donor_exon_reads'] / total_reads * yaxis_length
    df_sj_drawing['scaled_acceptor_exon_reads'] = df_sj_drawing['acceptor_exon_reads'] / total_reads * yaxis_length

    # draw splice junction reads. known donor: green; known acceptor: yellow; novel donor: purple; novel acceptor: orange
    for idx_label in df_sj_drawing[df_sj_drawing['is_known_donor'] == True].index:
        dwg.add(dwg.rect(insert=(df_sj_drawing.loc[idx_label, 'scaled_donor_exon_start'] + xoffset, exon_height + yoffset + 1), size=(df_sj_drawing.loc[idx_label, 'scaled_donor_exon_width'], df_sj_drawing.loc[idx_label, 'scaled_donor_exon_reads']), fill='lightgreen', fill_opacity=0.5))
    for idx_label in df_sj_drawing[df_sj_drawing['is_known_acceptor'] == True].index:
        dwg.add(dwg.rect(insert=(df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_start'] + xoffset, exon_height + yoffset + 1), size=(df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_width'], df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_reads']), fill='yellow', fill_opacity=0.5))
    for idx_label in df_sj_drawing[df_sj_drawing['is_known_donor'] == False].index:
        dwg.add(dwg.line(start=(df_sj_drawing.loc[idx_label, 'scaled_donor_exon_end'] + xoffset, exon_height + yoffset + 1), end=(df_sj_drawing.loc[idx_label, 'scaled_donor_exon_end'] + xoffset, df_sj_drawing.loc[idx_label, 'scaled_donor_exon_reads'] + exon_height + yoffset + 1), stroke='purple', stroke_width='0.25', stroke_opacity=1.0))
    for idx_label in df_sj_drawing[df_sj_drawing['is_known_acceptor'] == False].index:
        dwg.add(dwg.line(start=(df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_start'] + xoffset, exon_height + yoffset + 1), end=(df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_start'] + xoffset, df_sj_drawing.loc[idx_label, 'scaled_acceptor_exon_reads'] + exon_height + yoffset + 1), stroke='orange', stroke_width='0.25', stroke_opacity=1.0))
    # add title
    font_size = 5
    title_xstart = xaxis_length / 2 + xoffset - len(img_title) / 2 * font_size * 0.5
    dwg.add(dwg.text(img_title, insert=(title_xstart, yoffset / 2), fill='black', font_size=font_size))
    dwg.save()


os.chdir(path_splice_junction_folder)
list_sj_files = glob.glob(f'*{splice_junction_file_suffix}')

gtf_header_name = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'genomic_phase', 'info']
df_gtf = pd.read_csv(path_gtf_file, sep='\t', header=None, names=gtf_header_name, dtype={'chr': 'str', 'start': 'Int64', 'end': 'Int64', 'strand': 'str'}, low_memory=False)
df_gtf = df_gtf[df_gtf['feature_type'] == 'exon']
df_gtf = df_gtf[df_gtf['info'].str.contains(wanted)]
df_gtf = df_gtf.drop_duplicates(subset=['chr', 'start', 'end', 'strand'], keep='first')	 # drop duplicate exons except for the first occurrence
chr_name = list(set(df_gtf['chr']))[0]  # splice junction should be at the same chromosome, so do not need to iterate it
strand_name = list(set(df_gtf['strand']))[0]
if strand_name == '-':
    df_gtf['start'], df_gtf['end'] = -df_gtf['end'], -df_gtf['start']

df_all_sj = pd.DataFrame()  # save all splice junction data
splice_junction_header_name = [
    'chr', 'start', 'end', 'strand', 'intron_motif', 'is_annoted', 'number of uniquely mapping reads crossing the junction',
    'number of multi-mapping reads crossing the junction', 'maximum spliced alignment overhang']
for sj_file in list_sj_files:
    path_splice_junction_file = os.path.join(path_splice_junction_folder, sj_file)
    sj_filename_prefix = sj_file[:-len(splice_junction_file_suffix)].rstrip('.')
    path_output_svg_file = os.path.join(path_output_svg_folder, f'{sj_file}.svg')
    df_sj_data = pd.read_csv(path_splice_junction_file, sep='\t', header=None, names=splice_junction_header_name, dtype={'chr': 'str', 'start': 'Int64', 'end': 'Int64', 'strand': 'Int64'}, low_memory=False)
    df_sj_data['image_title'] = sj_filename_prefix
    df_sj_data['path_output_svg_file'] = path_output_svg_file
    df_all_sj = pd.concat([df_all_sj, df_sj_data], axis=0)
df_all_sj = df_all_sj.loc[df_all_sj['strand'] != 0]  # drop splice junctions on undefined strand
df_all_sj = df_all_sj.replace({'strand': {1: '+', 2: '-'}})
df_all_sj = df_all_sj[df_all_sj['strand'] == strand_name]
if strand_name == '-':
    df_all_sj['start'], df_all_sj['end'] = -df_all_sj['end'], -df_all_sj['start']
df_all_sj = df_all_sj[df_all_sj['number of uniquely mapping reads crossing the junction'] >= min_junction_reads_threshold]    # remove low junction reads
# select exons and splice junctions in the right ranges
df_all_sj = df_all_sj[(df_all_sj['chr'] == chr_name) & (df_all_sj['start'] > df_gtf['start'].min()) & (df_all_sj['end'] < df_gtf['end'].max())]
df_gtf = df_gtf[(df_gtf['end'] + 2 > df_all_sj['start'].min()) & (df_gtf['start'] < (df_all_sj['end'].max() + 2))]
df_all_sj = df_all_sj[(df_all_sj['chr'] == chr_name) & (df_all_sj['start'] > df_gtf['start'].min()) & (df_all_sj['end'] < df_gtf['end'].max())]
df_gtf = df_gtf.sort_values(['start', 'end'], ascending=True)
df_all_sj = df_all_sj.sort_values(['start', 'end'], ascending=True)

for path_output_svg_file in set(df_all_sj['path_output_svg_file'].tolist()):
    df_splice_junction = df_all_sj[df_all_sj['path_output_svg_file'] == path_output_svg_file].copy()
    img_title = df_splice_junction['image_title'].tolist()[0]
    df_splice_junction = df_splice_junction.drop(['image_title', 'path_output_svg_file'], axis=1)
    draw_sj_image(df_gtf, df_splice_junction, path_output_svg_file, img_title)
