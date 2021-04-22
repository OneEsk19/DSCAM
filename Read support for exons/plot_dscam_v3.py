
import svgwrite
import pandas as pd


min_junction_reads_threshold = 10

# get dscam exon coordinates from a GTF file
wanted = 'FBgn0033159'
filepath_gtf = '/home/user/dscam_plot/dscam_fruit_fly.gtf'  # I made me a smll gtf containing only Dscam, but would work on larger ones too (just slower)
filepath_junction = '/home/user/dscam_plot/final_script/NG-8304_Dscam1_lib80429_3962_GAATTCGTTATAGCCC_SJ.out.tab'
path_output_svg_file = '/home/user/dscam_plot/final_script/NG-8304_Dscam1_lib80429_3962_GAATTCGTTATAGCCC_SJ.svg'


gtf_header_name = ['chr', 'annotation_source', 'feature_type', 'start', 'end', 'score', 'strand', 'genomic_phase', 'info']
df_gtf = pd.read_csv(filepath_gtf, sep='\t', header=None, names=gtf_header_name, low_memory=False)
df_gtf = df_gtf[df_gtf['feature_type'] == 'exon']
df_gtf = df_gtf[df_gtf['info'].str.contains(wanted)]
df_gtf = df_gtf.drop_duplicates(subset=['start', 'end'], keep='first')	 # drop duplicate exons except for the first occurrence
chr_name = list(set(df_gtf['chr']))[0]  # splice junction should be at the same chromosome, so do not need to iterate it

# read the alignment information
splice_junction_header_name = [
    'chr', 'start', 'end', 'strand', 'intron_motif', 'is_annoted', 'number of uniquely mapping reads crossing the junction',
    'number of multi-mapping reads crossing the junction', 'maximum spliced alignment overhang']
df_splice_junction = pd.read_csv(filepath_junction, sep='\t', header=None, names=splice_junction_header_name, low_memory=False)
df_splice_junction = df_splice_junction[df_splice_junction['number of uniquely mapping reads crossing the junction'] > min_junction_reads_threshold]    # remove low junction reads
# select junctions in the right bondaries...
df_splice_junction = df_splice_junction[(df_splice_junction['chr'] == chr_name) & (df_splice_junction['start'] >= df_gtf['start'].min()) & (df_splice_junction['end'] <= df_gtf['end'].max())]

# select again
df_gtf = df_gtf[(df_gtf['end'] + 2 > df_splice_junction['start'].min()) & (df_gtf['start'] < df_splice_junction['end'].max() + 2)]
df_splice_junction = df_splice_junction[(df_splice_junction['chr'] == chr_name) & (df_splice_junction['start'] >= df_gtf['start'].min()) & (df_splice_junction['end'] <= df_gtf['end'].max())]
df_splice_junction = df_splice_junction.sort_values(['start', 'end'], ascending=True)

xoffset = 30
yoffset = 30
scale = 150
exon_height = 3

min_coordinate = df_gtf['start'].min()
max_coordinate = df_gtf['end'].max()

# draw the exons
dwg = svgwrite.Drawing(path_output_svg_file, profile='tiny')
for idx_label in df_gtf.index:
    start = (df_gtf.loc[idx_label, 'start'] - min_coordinate) / (max_coordinate - min_coordinate) * scale
    end = (df_gtf.loc[idx_label, 'end'] - min_coordinate) / (max_coordinate - min_coordinate) * scale
    width = (end - start)
    df_gtf.loc[idx_label, 'scaled_width'] = width
    dwg.add(dwg.rect(insert=(start + xoffset, yoffset), size=(width, exon_height), stroke='black', stroke_width=0.1, fill='lightblue', opacity=0.7))

df_splice_junction_start = pd.DataFrame()
df_splice_junction_end = pd.DataFrame()

for start_loc in set(df_splice_junction['start']):
    df_splice_junction_start.loc[start_loc, 'reads'] = df_splice_junction.loc[df_splice_junction['start']==start_loc, 'number of uniquely mapping reads crossing the junction'].sum()
for end_loc in set(df_splice_junction['end']):
    df_splice_junction_end.loc[end_loc, 'reads'] = df_splice_junction.loc[df_splice_junction['end']==end_loc, 'number of uniquely mapping reads crossing the junction'].sum()
df_splice_junction_start = df_splice_junction_start.sort_index()
df_splice_junction_end = df_splice_junction_end.sort_index()

max_reads_of_splice_junction = max(df_splice_junction_start['reads'].max(), df_splice_junction_end['reads'].max())
min_reads_of_splice_junction = min(df_splice_junction_start['reads'].min(), df_splice_junction_end['reads'].min())


def draw_junction_read(df_junction_reads, mode):
    if mode == 'start':
        for idx_label in df_junction_reads.index:
            start = (idx_label + 1 - min_coordinate) / (max_coordinate - min_coordinate) * scale
            width = df_gtf.loc[df_gtf['end'] == idx_label - 1, 'scaled_width'].tolist()[0]
            # height = (df_junction_reads.loc[idx_label, 'reads'] - min_reads_of_splice_junction) / (max_reads_of_splice_junction - min_reads_of_splice_junction) * scale
            height = df_junction_reads.loc[idx_label, 'reads'] / max_reads_of_splice_junction * scale
            dwg.add(dwg.rect(insert=(start - width + xoffset, exon_height + yoffset + 1), size=(width, height), stroke='black', stroke_width=0.1, stroke_opacity=0.5, fill='lightgreen', fill_opacity=0.7))
    elif mode == 'end':
        for idx_label in df_junction_reads.index:
            start = (idx_label + 1 - min_coordinate) / (max_coordinate - min_coordinate) * scale
            width = df_gtf.loc[df_gtf['start'] == idx_label + 1, 'scaled_width'].tolist()[0]
            # height = (df_junction_reads.loc[idx_label, 'reads'] - min_reads_of_splice_junction) / (max_reads_of_splice_junction - min_reads_of_splice_junction) * scale
            height = df_junction_reads.loc[idx_label, 'reads'] / max_reads_of_splice_junction * scale
            dwg.add(dwg.rect(insert=(start + xoffset, exon_height + yoffset + 1), size=(width, height), stroke='black', stroke_width=0.1, stroke_opacity=0.5, fill='yellow', fill_opacity=0.7))


draw_junction_read(df_splice_junction_start, 'start')
draw_junction_read(df_splice_junction_end, 'end')
dwg.save()
